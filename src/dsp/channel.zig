//! Generic channel impairment simulator: AWGN, wow/flutter (sinusoidal
//! time-warp), a dropout envelope, and gain — operates in place on a
//! caller-supplied slice with a caller-supplied `*std.Random`.
const std = @import("std");
const interp = @import("interp.zig");

const Sinc = interp.SincInterp(f32, 16, 64);

/// Add Gaussian noise so the result sits at `snr_db` relative to the
/// signal's measured RMS.
pub fn awgn(samples: []f32, snr_db: f32, rand: *std.Random) void {
    if (samples.len == 0) return;
    var power: f64 = 0;
    for (samples) |s| power += @as(f64, s) * @as(f64, s);
    power /= @floatFromInt(samples.len);
    const rms: f32 = @floatCast(@sqrt(power));
    const sigma = rms / std.math.pow(f32, 10.0, snr_db / 20.0);
    for (samples) |*s| s.* += sigma * rand.floatNorm(f32);
}

pub fn gain(samples: []f32, g: f32) void {
    for (samples) |*s| s.* *= g;
}

/// Parameters for `wowFlutter`.
pub const WowFlutter = struct {
    /// Slow pitch variation depth (fractional speed deviation) and rate (Hz).
    wow_depth: f32 = 0.002,
    wow_hz: f32 = 0.5,
    /// Fast mechanical variation depth (fractional) and rate (Hz).
    flutter_depth: f32 = 0.001,
    flutter_hz: f32 = 25.0,
    /// Constant playback speed error (fractional; +0.01 = 1% fast).
    speed_offset: f32 = 0.0,
    fs: f32 = 48_000.0,
};

/// Apply sinusoidal time-warp (wow/flutter, plus an optional constant speed
/// offset) via the windowed-sinc interpolator. `output.len` determines how
/// many warped samples are produced; for a nonzero `speed_offset` the
/// caller should size it as `floor(input.len / (1 + speed_offset))` to
/// avoid reading past the end of `input`. Overwrites `output` in place —
/// `input` and `output` must not alias.
pub fn wowFlutter(input: []const f32, output: []f32, p: WowFlutter) void {
    const fs64: f64 = @floatCast(p.fs);
    var read_pos: f64 = 0.0;
    for (0..output.len) |n| {
        const t: f64 = @as(f64, @floatFromInt(n)) / fs64;
        const speed: f64 = (1.0 + @as(f64, p.speed_offset)) +
            @as(f64, p.wow_depth) * @sin(2.0 * std.math.pi * @as(f64, p.wow_hz) * t) +
            @as(f64, p.flutter_depth) * @sin(2.0 * std.math.pi * @as(f64, p.flutter_hz) * t);
        output[n] = Sinc.sample(input, read_pos);
        read_pos += speed;
    }
}

/// Output length for `wowFlutter` given `in_len` and a constant speed
/// offset (0 for none): floor(in_len / (1 + speed_offset)).
pub fn wowFlutterOutLen(in_len: usize, speed_offset: f32) usize {
    if (speed_offset == 0) return in_len;
    return @intFromFloat(@floor(@as(f64, @floatFromInt(in_len)) / (1.0 + @as(f64, speed_offset))));
}

/// Parameters for `dropoutEnvelope`.
pub const Dropout = struct {
    /// Poisson event rate, events/sec; 0 disables.
    rate: f32 = 0.0,
    /// Duration range (uniform), seconds.
    min_s: f32 = 0.001,
    max_s: f32 = 0.010,
    /// Fade depth range (uniform), dB of attenuation.
    min_db: f32 = 10.0,
    max_db: f32 = 40.0,
    /// Raised-cosine edge ramp on each side of a dropout, seconds.
    edge_s: f32 = 0.0005,
    fs: f32 = 48_000.0,
};

/// Generate a per-sample gain envelope (1.0 = no dropout) into `env`
/// (length = the desired sample count). Multiply it into a signal to apply.
pub fn dropoutEnvelope(env: []f32, p: Dropout, rand: *std.Random) void {
    @memset(env, 1.0);
    if (p.rate <= 0) return;
    const n = env.len;
    const nf: f64 = @floatFromInt(n);
    var t: f64 = 0;
    while (true) {
        const u = @max(rand.float(f64), 1e-12);
        t += -@log(u) / @as(f64, p.rate);
        const start = t * @as(f64, p.fs);
        if (start >= nf) break;
        const dur_s = p.min_s + rand.float(f32) * (p.max_s - p.min_s);
        const depth_db = p.min_db + rand.float(f32) * (p.max_db - p.min_db);
        const g = std.math.pow(f32, 10.0, -depth_db / 20.0);
        const dur: f64 = @as(f64, dur_s) * @as(f64, p.fs);
        const edge: f32 = @max(p.edge_s * p.fs, 1.0);
        const lo: usize = @intFromFloat(@max(start, 0));
        const hi: usize = @min(n, @as(usize, @intFromFloat(start + dur)));
        for (lo..hi) |i| {
            const x: f32 = @floatCast(@as(f64, @floatFromInt(i)) - start);
            const from_end: f32 = @floatCast(start + dur - @as(f64, @floatFromInt(i)));
            var g_i: f32 = g;
            if (x < edge) {
                const c = 0.5 - 0.5 * @cos(std.math.pi * x / edge);
                g_i = 1.0 - (1.0 - g) * c;
            } else if (from_end < edge) {
                const c = 0.5 - 0.5 * @cos(std.math.pi * from_end / edge);
                g_i = 1.0 - (1.0 - g) * c;
            }
            env[i] = @min(env[i], g_i);
        }
    }
}

/// One example impairment bundle covering a typical lossy analog channel
/// (moderate wow/flutter, mild noise, occasional dropouts). Not tied to
/// any specific medium — a convenience starting point, not a preset table.
pub const example = struct {
    pub const snr_db: f32 = 45.0;
    pub const wow_flutter: WowFlutter = .{
        .wow_depth = 0.002,
        .wow_hz = 0.5,
        .flutter_depth = 0.001,
        .flutter_hz = 25.0,
    };
    pub const dropout: Dropout = .{
        .rate = 0.2,
        .min_s = 0.001,
        .max_s = 0.010,
        .min_db = 10.0,
        .max_db = 40.0,
    };
};

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test awgn {
    var prng = std.Random.DefaultPrng.init(1);
    var rand = prng.random();
    var sine: [8000]f32 = undefined;
    for (&sine, 0..) |*s, n| s.* = @sin(2.0 * std.math.pi * 1000.0 * @as(f32, @floatFromInt(n)) / 48_000.0);
    const orig = sine;
    awgn(&sine, 20.0, &rand);
    var noise_p: f64 = 0;
    var sig_p: f64 = 0;
    for (sine, orig) |a, b| {
        noise_p += @as(f64, a - b) * @as(f64, a - b);
        sig_p += @as(f64, b) * @as(f64, b);
    }
    const snr_measured = 10.0 * std.math.log10(sig_p / noise_p);
    try testing.expectApproxEqAbs(@as(f64, 20.0), snr_measured, 1.0);
}

test "awgn: no-op on empty slice" {
    var prng = std.Random.DefaultPrng.init(2);
    var rand = prng.random();
    var empty: [0]f32 = .{};
    awgn(&empty, 10.0, &rand);
}

test "gain: scales samples" {
    var s = [_]f32{ 1.0, -2.0, 0.5 };
    gain(&s, 2.0);
    try testing.expectEqualSlices(f32, &[_]f32{ 2.0, -4.0, 1.0 }, &s);
}

test "wowFlutter: zero depth and zero offset is (near) identity" {
    var input: [2000]f32 = undefined;
    for (&input, 0..) |*s, n| s.* = @sin(2.0 * std.math.pi * 1000.0 * @as(f32, @floatFromInt(n)) / 48_000.0);
    var out: [2000]f32 = undefined;
    wowFlutter(&input, &out, .{ .wow_depth = 0, .flutter_depth = 0, .speed_offset = 0 });
    for (input[32..1960], out[32..1960]) |a, b| try testing.expectApproxEqAbs(a, b, 1e-4);
}

test "wowFlutter: constant speed offset raises the zero-crossing rate" {
    const N = 4800;
    var input: [N]f32 = undefined;
    for (&input, 0..) |*s, n| s.* = @sin(2.0 * std.math.pi * 1000.0 * @as(f32, @floatFromInt(n)) / 48_000.0);
    const out_len = wowFlutterOutLen(N, 0.02);
    const out = try testing.allocator.alloc(f32, out_len);
    defer testing.allocator.free(out);
    wowFlutter(&input, out, .{ .wow_depth = 0, .flutter_depth = 0, .speed_offset = 0.02 });

    const H = struct {
        fn zcRate(s: []const f32) f64 {
            var c: usize = 0;
            for (1..s.len) |i| {
                if ((s[i - 1] < 0) != (s[i] < 0)) c += 1;
            }
            return @as(f64, @floatFromInt(c)) / @as(f64, @floatFromInt(s.len));
        }
    };
    const ratio = H.zcRate(out[0..4000]) / H.zcRate(input[0..4000]);
    try testing.expectApproxEqAbs(@as(f64, 1.02), ratio, 0.01);
    try testing.expectEqual(@as(usize, @intFromFloat(@floor(4800.0 / 1.02))), out_len);
}

test "dropoutEnvelope: unity baseline, events present, soft edges" {
    var prng = std.Random.DefaultPrng.init(7);
    var rand = prng.random();
    var env: [48_000]f32 = undefined;
    dropoutEnvelope(&env, .{ .rate = 10.0, .fs = 48_000.0 }, &rand);
    var min_v: f32 = 1.0;
    var max_step: f32 = 0;
    var dipped: usize = 0;
    for (env, 0..) |e, i| {
        try testing.expect(e <= 1.0 and e > 0.0);
        if (e < min_v) min_v = e;
        if (e < 0.99) dipped += 1;
        if (i > 0) max_step = @max(max_step, @abs(e - env[i - 1]));
    }
    try testing.expect(dipped > 0);
    try testing.expect(min_v < 0.4);
    try testing.expect(max_step <= 0.1);

    var prng0 = std.Random.DefaultPrng.init(7);
    var rand0 = prng0.random();
    var env0: [4_800]f32 = undefined;
    dropoutEnvelope(&env0, .{ .rate = 0, .fs = 48_000.0 }, &rand0);
    for (env0) |e| try testing.expectEqual(@as(f32, 1.0), e);
}

test "example preset builds and is usable" {
    var prng = std.Random.DefaultPrng.init(3);
    var rand = prng.random();
    var s = [_]f32{0.1} ** 100;
    awgn(&s, example.snr_db, &rand);
    var env: [100]f32 = undefined;
    dropoutEnvelope(&env, example.dropout, &rand);
    for (env) |e| try testing.expect(e > 0 and e <= 1.0);
}
