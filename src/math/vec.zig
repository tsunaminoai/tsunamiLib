const std = @import("std");

/// Free functions over native `@Vector` — arithmetic stays `a + b`, `a * s`
/// via `@splat`, and every op lowers to SIMD with no wrapper struct.
pub fn V(comptime n: comptime_int, comptime T: type) type {
    return @Vector(n, T);
}
pub const Vec2 = @Vector(2, f32);
pub const Vec3 = @Vector(3, f32);
pub const Vec4 = @Vector(4, f32);

fn Elem(comptime Vt: type) type {
    return @typeInfo(Vt).vector.child;
}
fn len_(comptime Vt: type) comptime_int {
    return @typeInfo(Vt).vector.len;
}

pub inline fn splat(comptime Vt: type, s: Elem(Vt)) Vt {
    return @splat(s);
}

pub inline fn scale(v: anytype, s: Elem(@TypeOf(v))) @TypeOf(v) {
    return v * @as(@TypeOf(v), @splat(s));
}

pub inline fn dot(a: anytype, b: @TypeOf(a)) Elem(@TypeOf(a)) {
    return @reduce(.Add, a * b);
}

pub inline fn lengthSq(v: anytype) Elem(@TypeOf(v)) {
    return dot(v, v);
}

pub inline fn length(v: anytype) Elem(@TypeOf(v)) {
    return @sqrt(dot(v, v));
}

pub inline fn distance(a: anytype, b: @TypeOf(a)) Elem(@TypeOf(a)) {
    return length(a - b);
}

/// Zero vector stays zero rather than producing NaN.
pub inline fn normalize(v: anytype) @TypeOf(v) {
    const l = length(v);
    return if (l == 0) v else scale(v, 1 / l);
}

pub inline fn lerp(a: anytype, b: @TypeOf(a), t: Elem(@TypeOf(a))) @TypeOf(a) {
    return a + scale(b - a, t);
}

pub inline fn clamp(v: anytype, lo: @TypeOf(v), hi: @TypeOf(v)) @TypeOf(v) {
    return @min(@max(v, lo), hi);
}

pub inline fn cross(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    comptime std.debug.assert(len_(@TypeOf(a)) == 3);
    const yzx = [3]i32{ 1, 2, 0 };
    const zxy = [3]i32{ 2, 0, 1 };
    return @shuffle(Elem(@TypeOf(a)), a, undefined, yzx) * @shuffle(Elem(@TypeOf(a)), b, undefined, zxy) -
        @shuffle(Elem(@TypeOf(a)), a, undefined, zxy) * @shuffle(Elem(@TypeOf(a)), b, undefined, yzx);
}

/// 2D perp-dot (z of the 3D cross).
pub inline fn cross2(a: anytype, b: @TypeOf(a)) Elem(@TypeOf(a)) {
    return a[0] * b[1] - a[1] * b[0];
}

/// Comptime swizzle: `swizzle(v, "zyx")`, `swizzle(v, "xxyy")`.
pub inline fn swizzle(v: anytype, comptime pattern: []const u8) @Vector(pattern.len, Elem(@TypeOf(v))) {
    const mask = comptime blk: {
        var m: [pattern.len]i32 = undefined;
        for (pattern, 0..) |c, i| m[i] = switch (c) {
            'x', 'r' => 0,
            'y', 'g' => 1,
            'z', 'b' => 2,
            'w', 'a' => 3,
            else => @compileError("bad swizzle component"),
        };
        for (m) |x| if (x >= len_(@TypeOf(v))) @compileError("swizzle out of range");
        break :blk m;
    };
    return @shuffle(Elem(@TypeOf(v)), v, undefined, mask);
}

pub inline fn extend(v: anytype, s: Elem(@TypeOf(v))) @Vector(len_(@TypeOf(v)) + 1, Elem(@TypeOf(v))) {
    const n = len_(@TypeOf(v));
    var out: @Vector(n + 1, Elem(@TypeOf(v))) = @splat(s);
    inline for (0..n) |i| out[i] = v[i];
    return out;
}

pub fn approxEq(a: anytype, b: @TypeOf(a), tol: Elem(@TypeOf(a))) bool {
    return @reduce(.And, @abs(a - b) <= @as(@TypeOf(a), @splat(tol)));
}

// ── Matrices ─────────────────────────────────────────────────────────────

/// Column-major n×n; `cols[j]` is column j so `mulVec` is n fused SIMD adds.
pub fn Mat(comptime n: comptime_int, comptime T: type) type {
    return struct {
        const Self = @This();
        pub const Col = @Vector(n, T);
        cols: [n]Col,

        pub const identity: Self = blk: {
            var m: Self = undefined;
            for (0..n) |j| {
                var c: Col = @splat(0);
                c[j] = 1;
                m.cols[j] = c;
            }
            break :blk m;
        };

        pub inline fn mulVec(m: Self, v: Col) Col {
            var acc: Col = m.cols[0] * @as(Col, @splat(v[0]));
            inline for (1..n) |j| acc += m.cols[j] * @as(Col, @splat(v[j]));
            return acc;
        }

        pub fn mul(a: Self, b: Self) Self {
            var r: Self = undefined;
            inline for (0..n) |j| r.cols[j] = a.mulVec(b.cols[j]);
            return r;
        }

        pub fn transpose(m: Self) Self {
            var r: Self = undefined;
            inline for (0..n) |i| inline for (0..n) |j| {
                r.cols[i][j] = m.cols[j][i];
            };
            return r;
        }

        pub fn at(m: Self, row: usize, col: usize) T {
            const c: [n]T = m.cols[col];
            return c[row];
        }
    };
}

pub const Mat3 = Mat(3, f32);
pub const Mat4 = Mat(4, f32);

pub fn translate(t: Vec3) Mat4 {
    var m = Mat4.identity;
    m.cols[3] = extend(t, 1);
    return m;
}

pub fn scaling(s: Vec3) Mat4 {
    var m = Mat4.identity;
    inline for (0..3) |i| m.cols[i][i] = s[i];
    return m;
}

/// Rodrigues rotation about a unit `axis` by `angle` radians (right-handed).
pub fn rotate(axis: Vec3, angle: f32) Mat4 {
    const a = normalize(axis);
    const c = @cos(angle);
    const s = @sin(angle);
    const t = 1 - c;
    return .{ .cols = .{
        .{ t * a[0] * a[0] + c, t * a[0] * a[1] + s * a[2], t * a[0] * a[2] - s * a[1], 0 },
        .{ t * a[0] * a[1] - s * a[2], t * a[1] * a[1] + c, t * a[1] * a[2] + s * a[0], 0 },
        .{ t * a[0] * a[2] + s * a[1], t * a[1] * a[2] - s * a[0], t * a[2] * a[2] + c, 0 },
        .{ 0, 0, 0, 1 },
    } };
}

/// Right-handed, clip z in [-1, 1] (OpenGL/raylib convention).
pub fn perspective(fovy: f32, aspect: f32, near: f32, far: f32) Mat4 {
    const f = 1 / @tan(fovy / 2);
    return .{ .cols = .{
        .{ f / aspect, 0, 0, 0 },
        .{ 0, f, 0, 0 },
        .{ 0, 0, (far + near) / (near - far), -1 },
        .{ 0, 0, 2 * far * near / (near - far), 0 },
    } };
}

pub fn lookAt(eye: Vec3, target: Vec3, up: Vec3) Mat4 {
    const f = normalize(target - eye);
    const s = normalize(cross(f, up));
    const u = cross(s, f);
    return .{ .cols = .{
        .{ s[0], u[0], -f[0], 0 },
        .{ s[1], u[1], -f[1], 0 },
        .{ s[2], u[2], -f[2], 0 },
        .{ -dot(s, eye), -dot(u, eye), dot(f, eye), 1 },
    } };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "basic vector ops" {
    const a: Vec3 = .{ 1, 2, 3 };
    const b: Vec3 = .{ 4, 5, 6 };
    try testing.expectEqual(@as(f32, 32), dot(a, b));
    try testing.expect(approxEq(cross(a, b), .{ -3, 6, -3 }, 1e-6));
    try testing.expectApproxEqAbs(@as(f32, 1), length(normalize(b)), 1e-6);
    try testing.expectEqual(@as(Vec3, @splat(0)), normalize(@as(Vec3, @splat(0))));
    try testing.expect(approxEq(lerp(a, b, 0.5), .{ 2.5, 3.5, 4.5 }, 1e-6));
    try testing.expectEqual(@as(f32, -3), cross2(Vec2{ 1, 2 }, Vec2{ 4, 5 }));
    try testing.expect(approxEq(clamp(b, @splat(0), @splat(5)), .{ 4, 5, 5 }, 0));
}

test "swizzle and extend" {
    const v: Vec4 = .{ 1, 2, 3, 4 };
    try testing.expectEqual(@Vector(3, f32){ 3, 2, 1 }, swizzle(v, "zyx"));
    try testing.expectEqual(@Vector(2, f32){ 4, 4 }, swizzle(v, "ww"));
    try testing.expectEqual(Vec4{ 1, 2, 3, 9 }, extend(Vec3{ 1, 2, 3 }, 9));
}

test "matrix identity, mul, transpose" {
    const t = translate(.{ 1, 2, 3 });
    const p = t.mulVec(.{ 1, 1, 1, 1 });
    try testing.expect(approxEq(p, .{ 2, 3, 4, 1 }, 1e-6));
    const ts = t.mul(scaling(.{ 2, 2, 2 }));
    try testing.expect(approxEq(ts.mulVec(.{ 1, 1, 1, 1 }), .{ 3, 4, 5, 1 }, 1e-6));
    try testing.expectEqual(@as(f32, 1), t.transpose().at(3, 0));
    const I = Mat4.identity;
    try testing.expectEqual(t, I.mul(t));
}

test "rotation and lookAt" {
    const r = rotate(.{ 0, 0, 1 }, std.math.pi / 2.0);
    try testing.expect(approxEq(r.mulVec(.{ 1, 0, 0, 1 }), .{ 0, 1, 0, 1 }, 1e-6));
    const view = lookAt(.{ 0, 0, 5 }, .{ 0, 0, 0 }, .{ 0, 1, 0 });
    try testing.expect(approxEq(view.mulVec(.{ 0, 0, 0, 1 }), .{ 0, 0, -5, 1 }, 1e-6));
    const proj = perspective(std.math.pi / 2.0, 1, 1, 100);
    const clip = proj.mulVec(.{ 0, 0, -1, 1 });
    try testing.expectApproxEqAbs(@as(f32, -1), clip[2] / clip[3], 1e-5);
}
