const std = @import("std");
const contract = @import("../contract.zig");

/// In-place iterative radix-2 FFT with comptime size. Twiddles and the
/// bit-reversal permutation are baked into rodata; stage loops are `inline for`
/// so every stride/half is a constant the optimizer can unroll against.
pub fn Fft(comptime T: type, comptime n: usize) type {
    comptime std.debug.assert(n >= 2 and std.math.isPowerOfTwo(n));
    return struct {
        pub const C = std.math.Complex(T);
        pub const len = n;
        const log2n = std.math.log2_int(usize, n);
        const Idx = std.math.IntFittingRange(0, n - 1);

        const twiddles: [n / 2]C = blk: {
            @setEvalBranchQuota(n * 64);
            var w: [n / 2]C = undefined;
            for (&w, 0..) |*t, k| {
                const a = -2.0 * std.math.pi * @as(f64, @floatFromInt(k)) / @as(f64, @floatFromInt(n));
                t.* = .{ .re = @floatCast(@cos(a)), .im = @floatCast(@sin(a)) };
            }
            break :blk w;
        };

        const bitrev: [n]Idx = blk: {
            @setEvalBranchQuota(n * 64);
            var r: [n]Idx = undefined;
            for (&r, 0..) |*v, i| v.* = @intCast(@bitReverse(@as(@Int(.unsigned, log2n), @intCast(i))));
            break :blk r;
        };

        pub fn forward(x: *[n]C) void {
            permute(x);
            inline for (0..log2n) |s| {
                const half = 1 << s;
                const stride = n >> (s + 1);
                var k: usize = 0;
                while (k < n) : (k += 2 * half) {
                    for (0..half) |j| {
                        const t = twiddles[j * stride].mul(x[k + j + half]);
                        const u = x[k + j];
                        x[k + j] = u.add(t);
                        x[k + j + half] = u.sub(t);
                    }
                }
            }
        }

        /// Scaled by 1/n so `inverse(forward(x)) == x`.
        pub fn inverse(x: *[n]C) void {
            for (x) |*v| v.im = -v.im;
            forward(x);
            const s: T = 1.0 / @as(T, n);
            for (x) |*v| v.* = .{ .re = v.re * s, .im = -v.im * s };
        }

        /// Real input, complex half-spectrum out (bins 0..n/2 inclusive).
        pub fn real(in: *const [n]T, scratch: *[n]C, out: []C) void {
            contract.require(out.len >= n / 2 + 1, "fft.real: out.len < n/2+1");
            for (scratch, in) |*c, r| c.* = .{ .re = r, .im = 0 };
            forward(scratch);
            @memcpy(out[0 .. n / 2 + 1], scratch[0 .. n / 2 + 1]);
        }

        inline fn permute(x: *[n]C) void {
            for (0..n) |i| {
                const j = bitrev[i];
                if (i < j) std.mem.swap(C, &x[i], &x[j]);
            }
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

fn naiveDft(comptime T: type, comptime n: usize, x: *const [n]std.math.Complex(T)) [n]std.math.Complex(T) {
    var out: [n]std.math.Complex(T) = undefined;
    for (&out, 0..) |*o, k| {
        o.* = .{ .re = 0, .im = 0 };
        for (x, 0..) |v, t| {
            const a = -2.0 * std.math.pi * @as(T, @floatFromInt(k * t)) / @as(T, n);
            o.* = o.add(v.mul(.{ .re = @cos(a), .im = @sin(a) }));
        }
    }
    return out;
}

test "matches naive DFT" {
    const F = Fft(f64, 64);
    var prng: std.Random.DefaultPrng = .init(1);
    const r = prng.random();
    var x: [64]F.C = undefined;
    for (&x) |*v| v.* = .{ .re = r.float(f64) - 0.5, .im = r.float(f64) - 0.5 };
    const want = naiveDft(f64, 64, &x);
    F.forward(&x);
    for (x, want) |g, w| {
        try testing.expectApproxEqAbs(w.re, g.re, 1e-9);
        try testing.expectApproxEqAbs(w.im, g.im, 1e-9);
    }
}

test "round trip and Parseval" {
    const F = Fft(f32, 1024);
    var prng: std.Random.DefaultPrng = .init(2);
    const r = prng.random();
    var x: [1024]F.C = undefined;
    for (&x) |*v| v.* = .{ .re = r.float(f32), .im = r.float(f32) };
    const orig = x;
    var e_t: f64 = 0;
    for (x) |v| e_t += v.re * v.re + v.im * v.im;
    F.forward(&x);
    var e_f: f64 = 0;
    for (x) |v| e_f += v.re * v.re + v.im * v.im;
    try testing.expectApproxEqRel(e_t, e_f / 1024.0, 1e-4);
    F.inverse(&x);
    for (x, orig) |g, w| try testing.expectApproxEqAbs(w.re, g.re, 1e-4);
}

test Fft {
    const F = Fft(f32, 256);
    var in: [256]f32 = undefined;
    for (&in, 0..) |*v, i| v.* = @cos(2.0 * std.math.pi * 10.0 * @as(f32, @floatFromInt(i)) / 256.0);
    var scratch: [256]F.C = undefined;
    var out: [129]F.C = undefined;
    F.real(&in, &scratch, &out);
    try testing.expectApproxEqAbs(@as(f32, 128), out[10].magnitude(), 1e-2);
    try testing.expectApproxEqAbs(@as(f32, 0), out[11].magnitude(), 1e-2);
}
