const std = @import("std");

pub fn gray(x: anytype) @TypeOf(x) {
    if (@bitSizeOf(@TypeOf(x)) < 2) return x;
    return x ^ (x >> 1);
}

pub fn ungray(x: anytype) @TypeOf(x) {
    const bits = @bitSizeOf(@TypeOf(x));
    var v = x;
    inline for (0..comptime std.math.log2_int_ceil(usize, @max(bits, 1))) |i| v ^= v >> (1 << i);
    return v;
}

/// Gray-coded PAM on one axis with levels ±1, ±3, … scaled by `scale`.
/// Code → level is a comptime table; level → code is O(1) slicing.
pub fn Pam(comptime T: type, comptime bits: u4, comptime scale: T) type {
    comptime std.debug.assert(bits >= 1);
    return struct {
        pub const Code = @Int(.unsigned, bits);
        pub const size = 1 << bits;

        /// levels[code]
        pub const levels: [size]T = blk: {
            var l: [size]T = undefined;
            for (0..size) |i| l[gray(@as(Code, @intCast(i)))] = (2 * @as(T, @floatFromInt(i)) - (size - 1)) * scale;
            break :blk l;
        };

        pub inline fn map(code: Code) T {
            return levels[code];
        }

        pub fn slice(x: T) Code {
            const pos = @round((x / scale + (size - 1)) / 2);
            const i: Code = @intFromFloat(std.math.clamp(pos, 0, size - 1));
            return gray(i);
        }

        /// Exact max-log LLRs, MSB first, positive = bit more likely 1.
        pub fn llr(x: T, noise_var: T, out: *[bits]T) void {
            var d0: [bits]T = @splat(std.math.floatMax(T));
            var d1: [bits]T = @splat(std.math.floatMax(T));
            inline for (0..size) |c| {
                const e = x - levels[c];
                const d = e * e;
                inline for (0..bits) |b| {
                    if ((c >> (bits - 1 - b)) & 1 == 1) d1[b] = @min(d1[b], d) else d0[b] = @min(d0[b], d);
                }
            }
            const inv = 1 / noise_var;
            for (out, d0, d1) |*o, a, b| o.* = (a - b) * inv;
        }
    };
}

/// Square QAM with unit average energy, built as Gray-PAM(I) × Gray-PAM(Q):
/// every grid neighbour differs in one bit and bit LLRs split exactly into
/// two independent PAM problems. Code layout: I bits high, Q bits low.
pub fn Qam(comptime T: type, comptime bits: u5) type {
    comptime std.debug.assert(bits >= 2 and bits % 2 == 0);
    const ab: u4 = bits / 2;
    const m: T = @floatFromInt(@as(u32, 1) << bits);
    const Axis = Pam(T, ab, 1 / @sqrt(2 * (m - 1) / 3));
    return struct {
        pub const C = std.math.Complex(T);
        pub const Code = @Int(.unsigned, bits);
        pub const size = 1 << bits;
        pub const axis = Axis;

        pub const points: [size]C = blk: {
            var p: [size]C = undefined;
            for (0..size) |c| p[c] = .{ .re = Axis.levels[c >> ab], .im = Axis.levels[c & (Axis.size - 1)] };
            break :blk p;
        };

        pub inline fn map(code: Code) C {
            return points[code];
        }

        pub fn slice(s: C) Code {
            return (@as(Code, Axis.slice(s.re)) << ab) | Axis.slice(s.im);
        }

        pub fn llr(s: C, noise_var: T, out: *[bits]T) void {
            Axis.llr(s.re, noise_var, out[0..ab]);
            Axis.llr(s.im, noise_var, out[ab..]);
        }

        /// Split a byte stream into codes; `out.len` must be ≥ bytes·8/bits.
        pub fn pack(bytes: []const u8, out: []Code) usize {
            var acc: u32 = 0;
            var have: u5 = 0;
            var n: usize = 0;
            for (bytes) |byte| {
                acc = (acc << 8) | byte;
                have += 8;
                while (have >= bits) : (n += 1) {
                    have -= bits;
                    out[n] = @truncate(acc >> have);
                }
                acc &= (@as(u32, 1) << have) - 1;
            }
            return n;
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "gray/ungray inverse" {
    for (0..256) |i| {
        const x: u8 = @intCast(i);
        try testing.expectEqual(x, ungray(gray(x)));
    }
    for (0..255) |i| {
        const a: u8 = @intCast(i);
        try testing.expectEqual(@as(u8, 1), @popCount(gray(a) ^ gray(a + 1)));
    }
}

test Qam {
    inline for (.{ 2, 4, 6, 8 }) |b| {
        const Q = Qam(f32, b);
        var e: f64 = 0;
        for (0..Q.size) |c| {
            const code: Q.Code = @intCast(c);
            try testing.expectEqual(code, Q.slice(Q.map(code)));
            const p = Q.map(code);
            e += p.re * p.re + p.im * p.im;
        }
        try testing.expectApproxEqAbs(1.0, e / Q.size, 1e-5);
    }
}

test "qam grid neighbours differ in one bit" {
    const Q = Qam(f64, 6);
    const step = Q.axis.levels[1] - Q.axis.levels[0];
    for (0..Q.size) |a| for (0..Q.size) |b| {
        const pa = Q.points[a];
        const pb = Q.points[b];
        const d = @abs(pa.re - pb.re) + @abs(pa.im - pb.im);
        if (@abs(d - @abs(step)) < 1e-9) try testing.expectEqual(@as(u7, 1), @popCount(a ^ b));
    };
}

test "llr signs match bits at the ideal point" {
    const Q = Qam(f32, 6);
    for (0..Q.size) |c| {
        var l: [6]f32 = undefined;
        Q.llr(Q.map(@intCast(c)), 0.1, &l);
        for (l, 0..) |v, b| {
            const bit = (c >> @intCast(5 - b)) & 1;
            try testing.expect((v > 0) == (bit == 1));
        }
    }
}

test "pack bytes into 6-bit codes" {
    const Q = Qam(f32, 6);
    var out: [4]Q.Code = undefined;
    try testing.expectEqual(@as(usize, 4), Q.pack(&.{ 0b111111_00, 0b0001_1010, 0b10_101010 }, &out));
    try testing.expectEqualSlices(Q.Code, &.{ 0b111111, 0b000001, 0b101010, 0b101010 }, &out);
}

test "slice clamps outliers" {
    const Q = Qam(f32, 4);
    try testing.expectEqual(Q.slice(Q.map(0b1010)), Q.slice(.{ .re = Q.map(0b1010).re * 50, .im = Q.map(0b1010).im * 50 }));
}
