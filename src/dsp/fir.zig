const std = @import("std");
const contract = @import("../contract.zig");
const window = @import("window.zig");

/// Comptime-designed windowed-sinc lowpass FIR, `taps` coefficients,
/// cutoff `fc` as a fraction of the sample rate (0 < fc < 0.5), normalised
/// to unity DC gain. `window` picks the taper: `.hann`, `.hamming`, `.blackman`.
pub fn lowpass(
    comptime T: type,
    comptime taps: usize,
    comptime fc: f64,
    comptime win: enum { hann, hamming, blackman },
) [taps]T {
    var h: [taps]T = undefined;
    const w = switch (win) {
        .hann => window.hann(f64, taps),
        .hamming => window.hamming(f64, taps),
        .blackman => window.blackman(f64, taps),
    };
    const centre: f64 = @as(f64, @floatFromInt(taps - 1)) / 2.0;
    var sum: f64 = 0;
    var proto: [taps]f64 = undefined;
    for (&proto, 0..) |*v, i| {
        const n = @as(f64, @floatFromInt(i)) - centre;
        const x = 2.0 * std.math.pi * fc * n;
        const sinc: f64 = if (@abs(n) < 1e-9) 1.0 else @sin(x) / x;
        v.* = 2.0 * fc * sinc * w[i];
        sum += v.*;
    }
    for (&h, proto) |*o, v| o.* = @floatCast(v / sum);
    return h;
}

/// Root-raised-cosine pulse shape, `sps` samples/symbol, `span` symbols
/// each side of centre (taps = 2·sps·span + 1), rolloff `beta` (0, 1].
/// Normalised to unity DC gain.
pub fn rootRaisedCosine(
    comptime T: type,
    comptime sps: usize,
    comptime span: usize,
    comptime beta: f64,
) [2 * sps * span + 1]T {
    const taps = 2 * sps * span + 1;
    var h: [taps]T = undefined;
    const centre: f64 = @as(f64, @floatFromInt(taps - 1)) / 2.0;
    var proto: [taps]f64 = undefined;
    var sum: f64 = 0;
    for (&proto, 0..) |*v, i| {
        const t = (@as(f64, @floatFromInt(i)) - centre) / @as(f64, @floatFromInt(sps));
        if (@abs(t) < 1e-9) {
            v.* = 1.0 - beta + 4.0 * beta / std.math.pi;
        } else if (@abs(@abs(4.0 * beta * t) - 1.0) < 1e-9) {
            const b = beta;
            v.* = (b / std.math.sqrt2) *
                ((1.0 + 2.0 / std.math.pi) * @sin(std.math.pi / (4.0 * b)) +
                    (1.0 - 2.0 / std.math.pi) * @cos(std.math.pi / (4.0 * b)));
        } else {
            const pit = std.math.pi * t;
            const num = @sin(pit * (1.0 - beta)) + 4.0 * beta * t * @cos(pit * (1.0 + beta));
            const denom = pit * (1.0 - (4.0 * beta * t) * (4.0 * beta * t));
            v.* = num / denom;
        }
        sum += v.*;
    }
    for (&h, proto) |*o, v| o.* = @floatCast(v / sum);
    return h;
}

/// Streaming FIR. The delay line is stored twice back-to-back so the newest
/// `taps` samples are always contiguous: one SIMD multiply + reduce per
/// sample, no modulo in the dot product.
pub fn Fir(comptime T: type, comptime taps: usize) type {
    const Vt = @Vector(taps, T);
    return struct {
        h_rev: Vt,
        line: [2 * taps]T = @splat(0),
        pos: usize = 0,

        const Self = @This();

        pub fn init(coeffs: [taps]T) Self {
            var r: [taps]T = undefined;
            for (&r, 0..) |*v, i| v.* = coeffs[taps - 1 - i];
            return .{ .h_rev = r };
        }

        pub fn process(self: *Self, x: T) T {
            self.line[self.pos] = x;
            self.line[self.pos + taps] = x;
            const recent: Vt = self.line[self.pos + 1 ..][0..taps].*;
            self.pos = if (self.pos + 1 == taps) 0 else self.pos + 1;
            return @reduce(.Add, self.h_rev * recent);
        }

        /// Filter a block into a caller-provided `out` slice (same length as `in`).
        pub fn processBlock(self: *Self, in: []const T, out: []T) void {
            contract.require(in.len == out.len, "Fir.processBlock: length mismatch");
            for (in, out) |x, *y| y.* = self.process(x);
        }

        pub fn reset(self: *Self) void {
            self.line = @splat(0);
            self.pos = 0;
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "lowpass: unity DC gain" {
    const h = lowpass(f32, 63, 0.1, .hamming);
    var sum: f32 = 0;
    for (h) |v| sum += v;
    try testing.expectApproxEqAbs(@as(f32, 1.0), sum, 1e-4);
}

test "lowpass: symmetric taps (linear phase)" {
    const h = lowpass(f64, 31, 0.15, .blackman);
    for (0..15) |i| try testing.expectApproxEqAbs(h[i], h[30 - i], 1e-12);
}

test "lowpass: attenuates a tone above cutoff, passes one below" {
    const N = 63;
    const h = lowpass(f32, N, 0.1, .hamming);
    var fir = Fir(f32, N).init(h);
    var out_hi: f32 = 0;
    var out_lo: f32 = 0;
    for (0..2000) |n| {
        const t: f32 = @floatFromInt(n);
        const hi = @sin(2.0 * std.math.pi * 0.3 * t);
        const yh = fir.process(hi);
        if (n > 1900) out_hi = @max(out_hi, @abs(yh));
    }
    fir.reset();
    for (0..2000) |n| {
        const t: f32 = @floatFromInt(n);
        const lo = @sin(2.0 * std.math.pi * 0.02 * t);
        const yl = fir.process(lo);
        if (n > 1900) out_lo = @max(out_lo, @abs(yl));
    }
    try testing.expect(out_hi < 0.05);
    try testing.expect(out_lo > 0.8);
}

test "rootRaisedCosine: unity DC gain and symmetry" {
    const h = rootRaisedCosine(f64, 4, 6, 0.25);
    var sum: f64 = 0;
    for (h) |v| sum += v;
    try testing.expectApproxEqAbs(@as(f64, 1.0), sum, 1e-6);
    const n = h.len;
    for (0..n / 2) |i| try testing.expectApproxEqAbs(h[i], h[n - 1 - i], 1e-9);
}

test "Fir: impulse response equals the taps" {
    const h = [_]f32{ 1, 2, 3, 4 };
    var fir = Fir(f32, 4).init(h);
    var out: [4]f32 = undefined;
    fir.processBlock(&[_]f32{ 1, 0, 0, 0 }, &out);
    try testing.expectEqualSlices(f32, &h, &out);
}

test "Fir: streaming process matches processBlock" {
    const h = [_]f32{ 0.25, 0.5, 0.25 };
    var a = Fir(f32, 3).init(h);
    var b = Fir(f32, 3).init(h);
    const in = [_]f32{ 1, 2, 3, 4, 5, 6 };
    var out_block: [in.len]f32 = undefined;
    b.processBlock(&in, &out_block);
    for (in, 0..) |x, i| {
        const y = a.process(x);
        try testing.expectApproxEqAbs(out_block[i], y, 1e-6);
    }
}
