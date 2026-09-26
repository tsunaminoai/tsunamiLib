const std = @import("std");
const contract = @import("../contract.zig");

/// Nearest entry parameter t >= 0 where a ray enters the axis-aligned box
/// `[bmin, bmax]`, or null on a miss (also null if the box is entirely
/// behind the ray). Origin inside the box reports t = 0. Slab test; `dir`
/// need not be normalized. `V` is `@Vector(n, T)`.
pub fn rayAabb(comptime V: type, orig: V, dir: V, bmin: V, bmax: V) ?std.meta.Child(V) {
    const T = std.meta.Child(V);
    const n = @typeInfo(V).vector.len;
    var t_near: T = 0.0;
    var t_far: T = std.math.inf(T);
    inline for (0..n) |axis| {
        const o = orig[axis];
        const d = dir[axis];
        if (d == 0) {
            if (o < bmin[axis] or o > bmax[axis]) return null;
        } else {
            const inv = 1.0 / d;
            var t0 = (bmin[axis] - o) * inv;
            var t1 = (bmax[axis] - o) * inv;
            if (t0 > t1) {
                const tmp = t0;
                t0 = t1;
                t1 = tmp;
            }
            t_near = @max(t_near, t0);
            t_far = @min(t_far, t1);
            if (t_near > t_far) return null;
        }
    }
    return t_near;
}

fn dot(comptime V: type, a: V, b: V) std.meta.Child(V) {
    return @reduce(.Add, a * b);
}

/// Squared distance from point `p` to segment `[a, b]`.
pub fn distSqPointSegment(comptime V: type, p: V, a: V, b: V) std.meta.Child(V) {
    const ab = b - a;
    const len2 = dot(V, ab, ab);
    const ap = p - a;
    if (len2 == 0) return dot(V, ap, ap);
    const t = std.math.clamp(dot(V, ap, ab) / len2, 0.0, 1.0);
    const closest = a + ab * @as(V, @splat(t));
    const d = p - closest;
    return dot(V, d, d);
}

/// Distance from point `p` to segment `[a, b]`.
pub fn distPointSegment(comptime V: type, p: V, a: V, b: V) std.meta.Child(V) {
    return @sqrt(distSqPointSegment(V, p, a, b));
}

/// Evaluate a cubic Bezier at parameter `t` in [0, 1]: `p0` start, `p1`/`p2`
/// control points, `p3` end.
pub fn bezierPoint(comptime V: type, t: std.meta.Child(V), p0: V, p1: V, p2: V, p3: V) V {
    const mt = 1.0 - t;
    const w0: std.meta.Child(V) = mt * mt * mt;
    const w1: std.meta.Child(V) = 3 * mt * mt * t;
    const w2: std.meta.Child(V) = 3 * mt * t * t;
    const w3: std.meta.Child(V) = t * t * t;
    return p0 * @as(V, @splat(w0)) + p1 * @as(V, @splat(w1)) +
        p2 * @as(V, @splat(w2)) + p3 * @as(V, @splat(w3));
}

/// Approximate distance from `p` to a cubic Bezier curve: minimum distance to
/// a fixed `steps`-segment polyline approximation. `steps` comptime so the
/// subdivision cost is chosen by the caller (12 for hit-testing, more for
/// precision).
pub fn distPointBezier(comptime V: type, comptime steps: usize, p: V, p0: V, p1: V, p2: V, p3: V) std.meta.Child(V) {
    comptime std.debug.assert(steps >= 1);
    const T = std.meta.Child(V);
    var best: T = std.math.inf(T);
    var prev = p0;
    for (1..steps + 1) |si| {
        const t: T = @as(T, @floatFromInt(si)) / @as(T, @floatFromInt(steps));
        const cur = bezierPoint(V, t, p0, p1, p2, p3);
        best = @min(best, distSqPointSegment(V, p, prev, cur));
        prev = cur;
    }
    return @sqrt(best);
}

/// Axis-aligned rectangle: `min`/`max` corners (min <= max on each axis).
pub fn Rect(comptime V: type) type {
    return struct {
        min: V,
        max: V,
    };
}

/// Whether `p` lies within `[r.min, r.max]` on every axis, inclusive.
pub fn rectContains(comptime V: type, r: Rect(V), p: V) bool {
    const n = @typeInfo(V).vector.len;
    inline for (0..n) |i| {
        if (p[i] < r.min[i] or p[i] > r.max[i]) return false;
    }
    return true;
}

/// Whether two axis-aligned rectangles overlap (touching edges count as
/// intersecting).
pub fn rectIntersect(comptime V: type, a: Rect(V), b: Rect(V)) bool {
    const n = @typeInfo(V).vector.len;
    inline for (0..n) |i| {
        if (a.max[i] < b.min[i] or b.max[i] < a.min[i]) return false;
    }
    return true;
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "rayAabb: hit, miss, origin inside, box behind" {
    const V = @Vector(3, f32);
    const bmin: V = .{ 1, -1, -1 };
    const bmax: V = .{ 3, 1, 1 };
    try testing.expectApproxEqAbs(@as(f32, 1.0), rayAabb(V, .{ 0, 0, 0 }, .{ 1, 0, 0 }, bmin, bmax).?, 1e-6);
    try testing.expect(rayAabb(V, .{ 0, 5, 0 }, .{ 1, 0, 0 }, bmin, bmax) == null);
    try testing.expectApproxEqAbs(@as(f32, 0.0), rayAabb(V, .{ 2, 0, 0 }, .{ 1, 0, 0 }, bmin, bmax).?, 1e-6);
    try testing.expect(rayAabb(V, .{ 5, 0, 0 }, .{ 1, 0, 0 }, bmin, bmax) == null);
    try testing.expect(rayAabb(V, .{ 0, 5, 0 }, .{ 0, 0, 1 }, bmin, bmax) == null);
}

test "distPointSegment: endpoints, midpoint, degenerate segment" {
    const V = @Vector(2, f32);
    const a: V = .{ 0, 0 };
    const b: V = .{ 10, 0 };
    try testing.expectApproxEqAbs(@as(f32, 0.0), distPointSegment(V, .{ 5, 0 }, a, b), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 3.0), distPointSegment(V, .{ 5, 3 }, a, b), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 5.0), distPointSegment(V, .{ -5, 0 }, a, b), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 5.0), distPointSegment(V, .{ 0, 5 }, a, a), 1e-6);
}

test "bezierPoint endpoints and midpoint of a symmetric curve" {
    const V = @Vector(2, f32);
    const p0: V = .{ 0, 0 };
    const p1: V = .{ 0, 10 };
    const p2: V = .{ 10, 10 };
    const p3: V = .{ 10, 0 };
    const start = bezierPoint(V, 0.0, p0, p1, p2, p3);
    const end = bezierPoint(V, 1.0, p0, p1, p2, p3);
    try testing.expectApproxEqAbs(@as(f32, 0), start[0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0), start[1], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 10), end[0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0), end[1], 1e-6);
    const mid = bezierPoint(V, 0.5, p0, p1, p2, p3);
    try testing.expectApproxEqAbs(@as(f32, 5), mid[0], 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 7.5), mid[1], 1e-6);
}

test "distPointBezier: zero on the curve, positive off it" {
    const V = @Vector(2, f32);
    const p0: V = .{ 0, 0 };
    const p1: V = .{ 0, 10 };
    const p2: V = .{ 10, 10 };
    const p3: V = .{ 10, 0 };
    const on_curve = bezierPoint(V, 0.3, p0, p1, p2, p3);
    try testing.expect(distPointBezier(V, 24, on_curve, p0, p1, p2, p3) < 0.05);
    try testing.expect(distPointBezier(V, 24, .{ -50, -50 }, p0, p1, p2, p3) > 40.0);
}

test "rectContains and rectIntersect" {
    const V = @Vector(2, f32);
    const R = Rect(V);
    const a: R = .{ .min = .{ 0, 0 }, .max = .{ 10, 10 } };
    const b: R = .{ .min = .{ 5, 5 }, .max = .{ 15, 15 } };
    const c: R = .{ .min = .{ 20, 20 }, .max = .{ 30, 30 } };
    try testing.expect(rectContains(V, a, .{ 5, 5 }));
    try testing.expect(rectContains(V, a, .{ 10, 10 }));
    try testing.expect(!rectContains(V, a, .{ 10.1, 5 }));
    try testing.expect(rectIntersect(V, a, b));
    try testing.expect(!rectIntersect(V, a, c));
}
