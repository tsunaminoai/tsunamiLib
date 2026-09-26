const std = @import("std");
const contract = @import("../contract.zig");

/// Modified Bessel function of the first kind, order 0 — power series
/// Σ ((x/2)^k / k!)², used by the Kaiser window.
pub fn besselI0(x: f64) f64 {
    const h = x / 2;
    var term: f64 = 1;
    var sum: f64 = 1;
    var k: f64 = 1;
    while (term > sum * 1e-17) : (k += 1) {
        const f = h / k;
        term *= f * f;
        sum += term;
    }
    return sum;
}

/// Hann window, n taps.
pub fn hann(comptime T: type, comptime n: usize) [n]T {
    var w: [n]T = undefined;
    for (&w, 0..) |*v, i| {
        const t: f64 = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n - 1));
        v.* = @floatCast(0.5 - 0.5 * @cos(2.0 * std.math.pi * t));
    }
    return w;
}

/// Hamming window, n taps.
pub fn hamming(comptime T: type, comptime n: usize) [n]T {
    var w: [n]T = undefined;
    for (&w, 0..) |*v, i| {
        const t: f64 = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n - 1));
        v.* = @floatCast(0.54 - 0.46 * @cos(2.0 * std.math.pi * t));
    }
    return w;
}

/// Blackman window, n taps.
pub fn blackman(comptime T: type, comptime n: usize) [n]T {
    var w: [n]T = undefined;
    for (&w, 0..) |*v, i| {
        const t: f64 = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n - 1));
        v.* = @floatCast(0.42 - 0.5 * @cos(2.0 * std.math.pi * t) + 0.08 * @cos(4.0 * std.math.pi * t));
    }
    return w;
}

/// Kaiser window, n taps, shape parameter beta (higher beta = more
/// sidelobe suppression, wider main lobe). beta ≈ 8 gives ≈ −81 dB sidelobes.
pub fn kaiser(comptime T: type, comptime n: usize, comptime beta: f64) [n]T {
    var w: [n]T = undefined;
    const norm = besselI0(beta);
    const centre: f64 = @as(f64, @floatFromInt(n - 1)) / 2.0;
    for (&w, 0..) |*v, i| {
        const r = (@as(f64, @floatFromInt(i)) - centre) / centre;
        const arg = beta * @sqrt(@max(0.0, 1.0 - r * r));
        v.* = @floatCast(besselI0(arg) / norm);
    }
    return w;
}

/// Multiply `window` into `samples` element-wise, in place.
pub fn apply(comptime T: type, window: []const T, samples: []T) void {
    contract.require(window.len == samples.len, "window.apply: length mismatch");
    for (samples, window) |*s, w| s.* *= w;
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "besselI0 matches known values" {
    try testing.expectApproxEqAbs(@as(f64, 1.0), besselI0(0.0), 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 1.2660658777520084), besselI0(1.0), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 427.56411572180474), besselI0(8.0), 1e-6);
}

test "hann window: zero endpoints, unity centre" {
    const w = hann(f64, 9);
    try testing.expectApproxEqAbs(@as(f64, 0.0), w[0], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 0.0), w[8], 1e-12);
    try testing.expectApproxEqAbs(@as(f64, 1.0), w[4], 1e-12);
}

test "hamming window: endpoints at 0.08, centre at 1.0" {
    const w = hamming(f64, 9);
    try testing.expectApproxEqAbs(@as(f64, 0.08), w[0], 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 1.0), w[4], 1e-9);
}

test "blackman window: endpoints near zero, centre at 1.0" {
    const w = blackman(f64, 9);
    try testing.expectApproxEqAbs(@as(f64, 0.0), w[0], 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 1.0), w[4], 1e-9);
}

test "kaiser window: symmetric, peak at centre" {
    const w = kaiser(f64, 15, 8.0);
    try testing.expectApproxEqAbs(@as(f64, 1.0), w[7], 1e-9);
    for (0..7) |i| try testing.expectApproxEqAbs(w[i], w[14 - i], 1e-9);
    try testing.expect(w[0] < w[7]);
}

test "apply scales samples in place" {
    const w = [_]f32{ 0.0, 1.0, 0.5 };
    var s = [_]f32{ 2.0, 3.0, 4.0 };
    apply(f32, &w, &s);
    try testing.expectEqualSlices(f32, &[_]f32{ 0.0, 3.0, 2.0 }, &s);
}
