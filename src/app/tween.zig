const std = @import("std");

/// Standard easing functions, `t` and result both in [0, 1].
pub const ease = struct {
    pub fn linear(t: f32) f32 {
        return t;
    }

    pub fn quadIn(t: f32) f32 {
        return t * t;
    }
    pub fn quadOut(t: f32) f32 {
        return t * (2.0 - t);
    }
    pub fn quadInOut(t: f32) f32 {
        if (t < 0.5) return 2.0 * t * t;
        const u = -2.0 * t + 2.0;
        return 1.0 - u * u / 2.0;
    }

    pub fn cubicIn(t: f32) f32 {
        return t * t * t;
    }
    pub fn cubicOut(t: f32) f32 {
        const u = 1.0 - t;
        return 1.0 - u * u * u;
    }
    pub fn cubicInOut(t: f32) f32 {
        if (t < 0.5) return 4.0 * t * t * t;
        const u = -2.0 * t + 2.0;
        return 1.0 - u * u * u / 2.0;
    }

    pub fn expoIn(t: f32) f32 {
        return if (t == 0) 0 else std.math.pow(f32, 2.0, 10.0 * t - 10.0);
    }
    pub fn expoOut(t: f32) f32 {
        return if (t == 1) 1 else 1.0 - std.math.pow(f32, 2.0, -10.0 * t);
    }
    pub fn expoInOut(t: f32) f32 {
        if (t == 0) return 0;
        if (t == 1) return 1;
        if (t < 0.5) return std.math.pow(f32, 2.0, 20.0 * t - 10.0) / 2.0;
        return (2.0 - std.math.pow(f32, 2.0, -20.0 * t + 10.0)) / 2.0;
    }

    // Overshoots past 1 before settling — "back" style, Penner's constants.
    pub fn backIn(t: f32) f32 {
        const c1: f32 = 1.70158;
        const c3: f32 = c1 + 1.0;
        return c3 * t * t * t - c1 * t * t;
    }
    pub fn backOut(t: f32) f32 {
        const c1: f32 = 1.70158;
        const c3: f32 = c1 + 1.0;
        const u = t - 1.0;
        return 1.0 + c3 * u * u * u + c1 * u * u;
    }
    pub fn backInOut(t: f32) f32 {
        const c2: f32 = 1.70158 * 1.525;
        if (t < 0.5) {
            const u = 2.0 * t;
            return (u * u * ((c2 + 1.0) * u - c2)) / 2.0;
        }
        const u = 2.0 * t - 2.0;
        return (u * u * ((c2 + 1.0) * u + c2) + 2.0) / 2.0;
    }

    pub fn elasticIn(t: f32) f32 {
        if (t == 0) return 0;
        if (t == 1) return 1;
        const c4: f32 = 2.0 * std.math.pi / 3.0;
        return -std.math.pow(f32, 2.0, 10.0 * t - 10.0) * @sin((t * 10.0 - 10.75) * c4);
    }
    pub fn elasticOut(t: f32) f32 {
        if (t == 0) return 0;
        if (t == 1) return 1;
        const c4: f32 = 2.0 * std.math.pi / 3.0;
        return std.math.pow(f32, 2.0, -10.0 * t) * @sin((t * 10.0 - 0.75) * c4) + 1.0;
    }
    pub fn elasticInOut(t: f32) f32 {
        if (t == 0) return 0;
        if (t == 1) return 1;
        const c5: f32 = 2.0 * std.math.pi / 4.5;
        if (t < 0.5) {
            return -(std.math.pow(f32, 2.0, 20.0 * t - 10.0) * @sin((20.0 * t - 11.125) * c5)) / 2.0;
        }
        return (std.math.pow(f32, 2.0, -20.0 * t + 10.0) * @sin((20.0 * t - 11.125) * c5)) / 2.0 + 1.0;
    }
};

/// A time-driven interpolation between `start` and `end` over `duration`
/// seconds, eased by comptime function `easeFn: fn(f32) f32`. `T` is any
/// value supporting `+`/`-`/scalar `*` — a float, or a `@Vector(n, f32)`.
pub fn Tween(comptime T: type, comptime easeFn: fn (f32) f32) type {
    return struct {
        start: T,
        end: T,
        duration: f32,
        /// Elapsed time, clamped to `[0, duration]`.
        elapsed: f32 = 0,

        const Self = @This();

        pub fn init(start: T, end: T, duration: f32) Self {
            return .{ .start = start, .end = end, .duration = duration };
        }

        /// Advance elapsed time by `dt`, clamped so it never runs past the end.
        pub fn update(self: *Self, dt: f32) void {
            self.elapsed = std.math.clamp(self.elapsed + dt, 0, self.duration);
        }

        pub fn done(self: Self) bool {
            return self.elapsed >= self.duration;
        }

        /// Fraction of duration elapsed, eased, in [0, 1].
        pub fn t(self: Self) f32 {
            const raw = if (self.duration <= 0) 1.0 else self.elapsed / self.duration;
            return easeFn(std.math.clamp(raw, 0, 1));
        }

        /// Current interpolated value.
        pub fn value(self: Self) T {
            const f = self.t();
            return lerp(self.start, self.end, f);
        }

        fn lerp(a: T, b: T, f: f32) T {
            return switch (@typeInfo(T)) {
                .float => a + (b - a) * @as(T, @floatCast(f)),
                .vector => |v| a + (b - a) * @as(T, @splat(@as(v.child, @floatCast(f)))),
                else => @compileError("Tween: unsupported value type " ++ @typeName(T)),
            };
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "linear tween reaches start, mid, and end" {
    var tw = Tween(f32, ease.linear).init(0, 10, 2.0);
    try testing.expect(!tw.done());
    try testing.expectApproxEqAbs(@as(f32, 0.0), tw.value(), 1e-6);
    tw.update(1.0);
    try testing.expectApproxEqAbs(@as(f32, 5.0), tw.value(), 1e-6);
    tw.update(1.0);
    try testing.expect(tw.done());
    try testing.expectApproxEqAbs(@as(f32, 10.0), tw.value(), 1e-6);
}

test "update clamps elapsed so it never overshoots duration" {
    var tw = Tween(f32, ease.linear).init(0, 4, 1.0);
    tw.update(100.0);
    try testing.expect(tw.done());
    try testing.expectApproxEqAbs(@as(f32, 4.0), tw.value(), 1e-6);
}

test "zero duration completes immediately" {
    var tw = Tween(f32, ease.quadOut).init(1, 2, 0.0);
    try testing.expect(tw.done());
    try testing.expectApproxEqAbs(@as(f32, 2.0), tw.value(), 1e-6);
}

test "vector value type interpolates componentwise" {
    const V = @Vector(2, f32);
    var tw = Tween(V, ease.linear).init(.{ 0, 0 }, .{ 10, -10 }, 2.0);
    tw.update(1.0);
    const v = tw.value();
    try testing.expectApproxEqAbs(@as(f32, 5.0), v[0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -5.0), v[1], 1e-6);
}

test "easing set: endpoints fixed at 0 and 1 for every curve" {
    const fns = [_]*const fn (f32) f32{
        ease.linear,    ease.quadIn,    ease.quadOut,    ease.quadInOut,
        ease.cubicIn,   ease.cubicOut,  ease.cubicInOut, ease.expoIn,
        ease.expoOut,   ease.expoInOut, ease.backIn,     ease.backOut,
        ease.backInOut, ease.elasticIn, ease.elasticOut, ease.elasticInOut,
    };
    for (fns) |f| {
        try testing.expectApproxEqAbs(@as(f32, 0.0), f(0.0), 1e-4);
        try testing.expectApproxEqAbs(@as(f32, 1.0), f(1.0), 1e-4);
    }
}

test "quadOut and quadIn are mirror images at the midpoint" {
    try testing.expectApproxEqAbs(@as(f32, 0.75), ease.quadOut(0.5), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.25), ease.quadIn(0.5), 1e-6);
}

test "backOut overshoots past 1 before settling" {
    var max: f32 = 0;
    var i: usize = 0;
    while (i <= 100) : (i += 1) {
        const tt: f32 = @as(f32, @floatFromInt(i)) / 100.0;
        max = @max(max, ease.backOut(tt));
    }
    try testing.expect(max > 1.0);
}
