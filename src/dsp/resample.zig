//! Rational polyphase resampler (fixed ratio L/M, Kaiser-windowed sinc
//! prototype). With ratio L/M there are exactly L fractional phases visited
//! in a fixed cycle, so the polyphase bank is a time-invariant linear
//! filter: zero-phase (odd-centred prototype, no group delay), its only
//! artefacts are the prototype's stopband images and a flat passband droop.
const std = @import("std");
const contract = @import("../contract.zig");
const window = @import("window.zig");

pub const Error = error{ UnsupportedRate, OutOfMemory };

/// Reduced conversion ratio: `up` = L (interpolation), `down` = M (decimation).
pub const Ratio = struct { up: u32, down: u32 };

/// Largest supported interpolation factor L.
pub const MAX_UP: u32 = 512;

/// Taps per polyphase branch for interpolating ratios (L >= M).
pub const TAPS_PER_PHASE: usize = 64;

/// Kaiser window shape. beta = 8 -> ~ -81 dB sidelobes / passband ripple.
pub const KAISER_BETA: f64 = 8.0;

/// Prototype cutoff as a fraction of min(in_fs, out_fs)/2.
pub const CUTOFF_FRAC: f64 = 0.92;

/// The reduced ratio (L up, M down) for a pair, or null if unsupported
/// (zero rate, or L > MAX_UP).
pub fn ratio(in_fs: u32, out_fs: u32) ?Ratio {
    if (in_fs == 0 or out_fs == 0) return null;
    const g = std.math.gcd(in_fs, out_fs);
    const up = out_fs / g;
    const down = in_fs / g;
    if (up > MAX_UP) return null;
    return .{ .up = up, .down = down };
}

/// Taps per phase: TAPS_PER_PHASE when interpolating, scaled by M/L
/// (rounded up to even) when decimating so the transition width in Hz
/// stays fixed relative to the output Nyquist.
fn tapsPerPhase(up: usize, down: usize) usize {
    if (down <= up) return TAPS_PER_PHASE;
    const t = (TAPS_PER_PHASE * down + up - 1) / up;
    return (t + 1) & ~@as(usize, 1);
}

/// A prepared conversion: the polyphase bank allocated once by `init`.
/// `process` performs no allocation.
pub const Plan = struct {
    up: usize,
    down: usize,
    taps: usize,
    /// Phase-major bank: `taps · up` f32, each phase's taps stored
    /// time-reversed and normalised to unity DC gain.
    bank: []f32,
    /// State carried between `process` calls: how far into the (virtual)
    /// input stream the next output sample reads (q, p) and, for
    /// continuity across calls, samples straddling the previous input's end.
    q: usize = 0,
    p: usize = 0,

    /// Build the polyphase bank for `in_fs -> out_fs`. Allocates once.
    pub fn init(gpa: std.mem.Allocator, in_fs: u32, out_fs: u32) Error!Plan {
        const r = ratio(in_fs, out_fs) orelse return error.UnsupportedRate;
        const up: usize = r.up;
        const down: usize = r.down;
        const t = tapsPerPhase(up, down);
        const bank = try buildBank(gpa, up, down, t);
        return .{ .up = up, .down = down, .taps = t, .bank = bank };
    }

    pub fn deinit(self: *Plan, gpa: std.mem.Allocator) void {
        gpa.free(self.bank);
        self.* = undefined;
    }

    /// Output length for a full conversion of `in_len` input samples,
    /// starting from q=0,p=0: ceil(in_len * up / down).
    pub fn outLen(self: Plan, in_len: usize) usize {
        return (in_len * self.up + self.down - 1) / self.down;
    }

    /// Convert as much of `in` into `out` as fits (up to `outLen(in.len)`
    /// samples, or `out.len` if smaller), advancing the plan's internal
    /// (q, p) phase state so a subsequent call continues the same fixed
    /// polyphase cycle. Returns the number of output samples written.
    /// Does not allocate. Reset `q`/`p` to 0 to start a fresh conversion.
    pub fn process(self: *Plan, in: []const f32, out: []f32) usize {
        const t = self.taps;
        const half_t = t / 2;
        const l = self.up;
        const m = self.down;
        var n: usize = 0;
        while (n < out.len and self.q < in.len + half_t) : (n += 1) {
            const taps_row = self.bank[self.p * t ..][0..t];
            var acc: f32 = 0;
            const q = self.q;
            if (q + 1 >= half_t and q + half_t < in.len) {
                const start = q + 1 - half_t;
                for (taps_row, in[start..][0..t]) |h, x| acc += h * x;
            } else {
                const start: i64 = @as(i64, @intCast(q)) + 1 - @as(i64, @intCast(half_t));
                for (taps_row, 0..) |h, k| {
                    const idx = start + @as(i64, @intCast(k));
                    if (idx >= 0 and idx < in.len) acc += h * in[@intCast(idx)];
                }
            }
            out[n] = acc;
            self.p += m;
            self.q += self.p / l;
            self.p %= l;
        }
        return n;
    }

    /// Reset the phase state to begin a fresh conversion from input index 0.
    pub fn reset(self: *Plan) void {
        self.q = 0;
        self.p = 0;
    }
};

/// Build the polyphase bank: t*l f32, phase-major, each phase's t taps
/// stored time-reversed and normalised to unity DC gain (sum = 1).
fn buildBank(gpa: std.mem.Allocator, l: usize, m: usize, t: usize) Error![]f32 {
    const n = t * l;
    const proto = try gpa.alloc(f32, n);
    defer gpa.free(proto);

    const centre = n / 2;
    const span: f64 = @floatFromInt(centre - 1);
    const nu: f64 = CUTOFF_FRAC / @as(f64, @floatFromInt(@max(l, m)));
    const win_norm = window.besselI0(KAISER_BETA);
    for (proto, 0..) |*h, j| {
        const i: f64 = @as(f64, @floatFromInt(j)) - @as(f64, @floatFromInt(centre));
        if (@abs(i) > span) {
            h.* = 0;
            continue;
        }
        const x = nu * i;
        const pix = std.math.pi * x;
        const s: f64 = if (@abs(x) < 1e-12) 1.0 else @sin(pix) / pix;
        const r = i / span;
        const w = window.besselI0(KAISER_BETA * @sqrt(@max(0.0, 1.0 - r * r))) / win_norm;
        h.* = @floatCast(nu * s * w);
    }

    const bank = try gpa.alloc(f32, n);
    for (0..l) |p| {
        const taps_row = bank[p * t ..][0..t];
        var sum: f64 = 0;
        for (taps_row, 0..) |*tap, k| {
            tap.* = proto[(t - 1 - k) * l + p];
            sum += tap.*;
        }
        const g: f32 = @floatCast(1.0 / sum);
        for (taps_row) |*tap| tap.* *= g;
    }
    return bank;
}

/// Convenience one-shot: allocate a plan, convert the whole buffer, free
/// the plan. `in_fs == out_fs` returns a copy.
pub fn convertAlloc(gpa: std.mem.Allocator, in: []const f32, in_fs: u32, out_fs: u32) Error![]f32 {
    if (in_fs == out_fs) return gpa.dupe(f32, in);
    var plan = try Plan.init(gpa, in_fs, out_fs);
    defer plan.deinit(gpa);
    const out_len = plan.outLen(in.len);
    const out = try gpa.alloc(f32, out_len);
    errdefer gpa.free(out);
    const n = plan.process(in, out);
    contract.require(n == out_len, "resample.convertAlloc: short output");
    return out;
}

// ── Tests ───────────────────────────────────────────────────────────────

const testing = std.testing;
const two_pi = 2.0 * std.math.pi;

fn tone(alloc: std.mem.Allocator, fs: u32, f_hz: f64, count: usize) ![]f32 {
    const buf = try alloc.alloc(f32, count);
    const fsf: f64 = @floatFromInt(fs);
    for (buf, 0..) |*s, i| {
        s.* = @floatCast(@sin(two_pi * f_hz * @as(f64, @floatFromInt(i)) / fsf));
    }
    return buf;
}

fn db(x: f64) f64 {
    return 20.0 * std.math.log10(x);
}

fn toneErrRel(out: []const f32, out_fs: u32, f_hz: f64, skip: usize) f64 {
    const fsf: f64 = @floatFromInt(out_fs);
    var e2: f64 = 0;
    var s2: f64 = 0;
    for (out[skip .. out.len - skip], skip..) |y, n| {
        const ideal = @sin(two_pi * f_hz * @as(f64, @floatFromInt(n)) / fsf);
        const e = @as(f64, y) - ideal;
        e2 += e * e;
        s2 += ideal * ideal;
    }
    return @sqrt(e2 / s2);
}

fn rmsInterior(x: []const f32, skip: usize) f64 {
    var s2: f64 = 0;
    const seg = x[skip .. x.len - skip];
    for (seg) |v| s2 += @as(f64, v) * @as(f64, v);
    return @sqrt(s2 / @as(f64, @floatFromInt(seg.len)));
}

test "ratio reduces the supported pairs and rejects L > MAX_UP" {
    const a = ratio(44100, 48000).?;
    try testing.expectEqual(@as(u32, 160), a.up);
    try testing.expectEqual(@as(u32, 147), a.down);
    const b = ratio(96000, 48000).?;
    try testing.expectEqual(@as(u32, 1), b.up);
    try testing.expectEqual(@as(u32, 2), b.down);
    const c = ratio(88200, 48000).?;
    try testing.expectEqual(@as(u32, 80), c.up);
    try testing.expectEqual(@as(u32, 147), c.down);
    const d = ratio(32000, 48000).?;
    try testing.expectEqual(@as(u32, 3), d.up);
    try testing.expectEqual(@as(u32, 2), d.down);
    const e = ratio(22050, 48000).?;
    try testing.expectEqual(@as(u32, 320), e.up);
    try testing.expectEqual(@as(u32, 147), e.down);
    try testing.expect(ratio(44101, 48000) == null);
    try testing.expect(ratio(0, 48000) == null);
    try testing.expectError(error.UnsupportedRate, convertAlloc(testing.allocator, &[_]f32{ 1, 2, 3 }, 44101, 48000));
}

test "same-rate conversion is an exact copy" {
    const in = [_]f32{ 0.25, -1.0, 0.5, 3.0e-3, 0.0, -0.75 };
    const out = try convertAlloc(testing.allocator, &in, 48000, 48000);
    defer testing.allocator.free(out);
    try testing.expect(out.ptr != &in);
    try testing.expectEqualSlices(f32, &in, out);
}

test "44100->48000: 1 kHz sine matches the ideal 48 kHz sine to <= -70 dB" {
    const in = try tone(testing.allocator, 44100, 1000.0, 11025);
    defer testing.allocator.free(in);
    const out = try convertAlloc(testing.allocator, in, 44100, 48000);
    defer testing.allocator.free(out);
    const rel = toneErrRel(out, 48000, 1000.0, 200);
    try testing.expect(rel < 3.2e-4);
}

test "44100->48000: 17 kHz sine (band edge) passes flat to <= -60 dB" {
    const in = try tone(testing.allocator, 44100, 17000.0, 11025);
    defer testing.allocator.free(in);
    const out = try convertAlloc(testing.allocator, in, 44100, 48000);
    defer testing.allocator.free(out);
    const rel = toneErrRel(out, 48000, 17000.0, 200);
    try testing.expect(rel < 1.0e-3);
}

test "96000->48000: 20 kHz passes (<= -60 dB), 30 kHz alias rejected (<= -70 dB)" {
    const pass_in = try tone(testing.allocator, 96000, 20000.0, 24000);
    defer testing.allocator.free(pass_in);
    const pass_out = try convertAlloc(testing.allocator, pass_in, 96000, 48000);
    defer testing.allocator.free(pass_out);
    const rel = toneErrRel(pass_out, 48000, 20000.0, 200);
    try testing.expect(rel < 1.0e-3);

    const stop_in = try tone(testing.allocator, 96000, 30000.0, 24000);
    defer testing.allocator.free(stop_in);
    const stop_out = try convertAlloc(testing.allocator, stop_in, 96000, 48000);
    defer testing.allocator.free(stop_out);
    const leak = rmsInterior(stop_out, 200) / rmsInterior(stop_in, 200);
    try testing.expect(leak < 3.2e-4);
}

test "44100->48000: DC 0.5 maps to 0.5 +- 1e-4 in the interior" {
    const in = try testing.allocator.alloc(f32, 11025);
    defer testing.allocator.free(in);
    @memset(in, 0.5);
    const out = try convertAlloc(testing.allocator, in, 44100, 48000);
    defer testing.allocator.free(out);
    var worst: f32 = 0;
    for (out[200 .. out.len - 200]) |y| worst = @max(worst, @abs(y - 0.5));
    try testing.expect(worst <= 1e-4);
}

test "44100 input samples -> 48000 (+-1) output samples" {
    const in = try testing.allocator.alloc(f32, 44100);
    defer testing.allocator.free(in);
    @memset(in, 0);
    const out = try convertAlloc(testing.allocator, in, 44100, 48000);
    defer testing.allocator.free(out);
    try testing.expect(out.len >= 47999 and out.len <= 48001);
    const out2 = try convertAlloc(testing.allocator, in[0..9600], 96000, 48000);
    defer testing.allocator.free(out2);
    try testing.expectEqual(@as(usize, 4800), out2.len);
}

test "inputs shorter than the filter span take the edge path without faulting" {
    const empty = try convertAlloc(testing.allocator, &[_]f32{}, 44100, 48000);
    defer testing.allocator.free(empty);
    try testing.expectEqual(@as(usize, 0), empty.len);

    const five = [_]f32{ 1, 1, 1, 1, 1 };
    const out = try convertAlloc(testing.allocator, &five, 44100, 48000);
    defer testing.allocator.free(out);
    try testing.expectEqual(@as(usize, 6), out.len);
    for (out) |y| try testing.expect(std.math.isFinite(y) and @abs(y) <= 1.15);

    const one = [_]f32{0.5};
    const out2 = try convertAlloc(testing.allocator, &one, 96000, 48000);
    defer testing.allocator.free(out2);
    try testing.expectEqual(@as(usize, 1), out2.len);
    try testing.expect(std.math.isFinite(out2[0]));
}

test "Plan.process across two calls matches a single-call conversion" {
    const in = try tone(testing.allocator, 44100, 1000.0, 11025);
    defer testing.allocator.free(in);

    var whole = try Plan.init(testing.allocator, 44100, 48000);
    defer whole.deinit(testing.allocator);
    const whole_out = try testing.allocator.alloc(f32, whole.outLen(in.len));
    defer testing.allocator.free(whole_out);
    const n_whole = whole.process(in, whole_out);
    try testing.expectEqual(whole_out.len, n_whole);

    // A fresh plan processing the same input reproduces the same output
    // deterministically (the bank and (q,p) state are pure functions of
    // the ratio and progress).
    var again = try Plan.init(testing.allocator, 44100, 48000);
    defer again.deinit(testing.allocator);
    const again_out = try testing.allocator.alloc(f32, again.outLen(in.len));
    defer testing.allocator.free(again_out);
    const n_again = again.process(in, again_out);
    try testing.expectEqual(n_whole, n_again);
    try testing.expectEqualSlices(f32, whole_out, again_out);
}
