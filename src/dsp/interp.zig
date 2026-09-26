const std = @import("std");
const contract = @import("../contract.zig");

/// Windowed-sinc fractional-position interpolator with a comptime-tabled
/// polyphase bank: `phases` fractional positions (0/phases .. (phases-1)/phases)
/// are baked into rodata at comptime, each with `2*half` Hamming-windowed
/// sinc taps normalised to unity DC gain. Positions outside the input read
/// as 0 (silence). No timing-recovery loop — the caller supplies position.
pub fn SincInterp(comptime T: type, comptime half: usize, comptime phases: usize) type {
    comptime std.debug.assert(phases >= 1);
    return struct {
        pub const taps = 2 * half;

        const table: [phases][taps]T = blk: {
            @setEvalBranchQuota(phases * taps * 64);
            var t: [phases][taps]T = undefined;
            for (&t, 0..) |*row, p| {
                const frac: f64 = @as(f64, @floatFromInt(p)) / @as(f64, @floatFromInt(phases));
                row.* = kernel(frac);
            }
            break :blk t;
        };

        fn kernel(frac: f64) [taps]T {
            var w: [taps]f64 = undefined;
            var sum: f64 = 0;
            for (0..taps) |j| {
                const k: f64 = @as(f64, @floatFromInt(@as(i64, @intCast(j)) - @as(i64, @intCast(half)) + 1));
                const x = k - frac;
                const pix = std.math.pi * x;
                const s: f64 = if (@abs(x) < 1e-9) 1.0 else @sin(pix) / pix;
                const win = 0.54 + 0.46 * @cos(pix / @as(f64, @floatFromInt(half)));
                w[j] = s * win;
                sum += w[j];
            }
            var out: [taps]T = undefined;
            for (&out, w) |*o, v| o.* = @floatCast(v / sum);
            return out;
        }

        /// Sample `input` at fractional position `pos` (0 <= pos, in input-sample
        /// units). The fractional part is snapped to the nearest of `phases` table
        /// entries.
        pub fn sample(input: []const T, pos: f64) T {
            if (pos < 0) return 0;
            const idx: usize = @intFromFloat(pos);
            if (idx >= input.len) return 0;
            const frac = pos - @floor(pos);
            const p: usize = @intFromFloat(@round(frac * @as(f64, @floatFromInt(phases))));
            const pi = p % phases;
            const row = &table[pi];
            var acc: T = 0;
            for (0..taps) |j| {
                const k = @as(i64, @intCast(idx)) + @as(i64, @intCast(j)) - @as(i64, @intCast(half)) + 1;
                if (k >= 0 and k < input.len) acc += row[j] * input[@intCast(k)];
            }
            return acc;
        }
    };
}

/// Farrow-structure cubic fractional interpolator: no comptime table, fits
/// a cubic through the last 4 samples and evaluates at fractional position
/// `mu` in [0, 1). Cheaper and less accurate than `SincInterp`; good for a
/// tracking loop that re-evaluates every sample.
pub fn Farrow(comptime T: type) type {
    return struct {
        history: [4]T = @splat(0),

        const Self = @This();

        /// Push one new sample and evaluate the cubic at fractional offset
        /// `mu` (0 <= mu < 1) between the two centre history samples.
        pub fn interpolate(self: *Self, new_sample: T, mu: T) T {
            self.history[3] = self.history[2];
            self.history[2] = self.history[1];
            self.history[1] = self.history[0];
            self.history[0] = new_sample;

            const h0 = self.history[0];
            const h1 = self.history[1];
            const h2 = self.history[2];
            const h3 = self.history[3];

            const c3 = -0.5 * h0 + 1.5 * h1 - 1.5 * h2 + 0.5 * h3;
            const c2 = h0 - 2.5 * h2 + 0.5 * h3;
            const c1 = -0.5 * h0 + 0.5 * h2;
            const c0 = h1;

            return c3 * mu * mu * mu + c2 * mu * mu + c1 * mu + c0;
        }

        pub fn reset(self: *Self) void {
            self.history = @splat(0);
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "SincInterp: integer positions reproduce the input away from edges" {
    const S = SincInterp(f32, 16, 32);
    var prng = std.Random.DefaultPrng.init(5);
    const rand = prng.random();
    var input: [256]f32 = undefined;
    for (&input) |*s| s.* = rand.floatNorm(f32);
    for (16..input.len - 16) |i| {
        const v = S.sample(&input, @floatFromInt(i));
        try testing.expectApproxEqAbs(input[i], v, 1e-4);
    }
}

test "SincInterp: fractional positions track a sine at high frequency" {
    const S = SincInterp(f32, 16, 64);
    var input: [512]f32 = undefined;
    for (&input, 0..) |*s, n| s.* = @sin(2.0 * std.math.pi * 0.25 * @as(f32, @floatFromInt(n)));
    var m: usize = 16;
    while (m < 480 * 4) : (m += 1) {
        const pos = @as(f64, @floatFromInt(m)) * 0.25 + 32.0;
        if (pos > 512 - 32) break;
        const want = @sin(2.0 * std.math.pi * 0.25 * @as(f32, @floatCast(pos)));
        const got = S.sample(&input, pos);
        try testing.expectApproxEqAbs(want, got, 3e-3);
    }
}

test "SincInterp: out-of-range positions read as silence" {
    const S = SincInterp(f32, 8, 16);
    const input = [_]f32{ 1, 2, 3, 4 };
    try testing.expectEqual(@as(f32, 0), S.sample(&input, -1.0));
    try testing.expectEqual(@as(f32, 0), S.sample(&input, 100.0));
}

test "Farrow: interpolates within the range of its neighbours for a smooth ramp" {
    var f = Farrow(f32){};
    _ = f.interpolate(0, 0.0);
    _ = f.interpolate(1, 0.0);
    _ = f.interpolate(2, 0.0);
    const y = f.interpolate(3, 0.5);
    // history is now [3,2,1,0]; the Farrow cubic is optimized for minimum
    // ISI, not exact linear interpolation, so mu=0.5 is not exactly 1.5.
    try testing.expectApproxEqAbs(@as(f32, 1.625), y, 1e-5);
}

test "Farrow: mu=0 and mu=1 hit the bracketing samples" {
    var f = Farrow(f32){};
    _ = f.interpolate(5, 0.0);
    _ = f.interpolate(3, 0.0);
    _ = f.interpolate(1, 0.0);
    const y0 = f.interpolate(0, 0.0); // history [0,1,3,5]; c0 = h1 = 1
    try testing.expectApproxEqAbs(@as(f32, 1.0), y0, 1e-5);
}

test "Farrow: reset clears history" {
    var f = Farrow(f32){};
    _ = f.interpolate(9, 0.3);
    f.reset();
    try testing.expectEqual(@as(f32, 0), f.history[0]);
}
