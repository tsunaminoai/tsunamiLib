//! Tour of math.vec (cross/swizzle/Mat4), math.geom (rayAabb, bezier
//! distance), math.geo (great-circle distance/bearing), math.grid
//! (neighbor iteration) and math.astro (GMST at a given date).
const std = @import("std");
const ts = @import("tsunami");
const vec = ts.math.vec;
const geom = ts.math.geom;
const geo = ts.math.geo;
const grid = ts.math.grid;
const astro = ts.math.astro;

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    // ── vec: cross, swizzle, Mat4 project ──
    const a: vec.Vec3 = .{ 1, 0, 0 };
    const b: vec.Vec3 = .{ 0, 1, 0 };
    const up = vec.cross(a, b);
    const swz = vec.swizzle(vec.Vec4{ 1, 2, 3, 4 }, "wzyx");

    const view = vec.lookAt(.{ 0, 0, 5 }, .{ 0, 0, 0 }, .{ 0, 1, 0 });
    const proj = vec.perspective(std.math.pi / 2.0, 1.0, 0.1, 100.0);
    const clip = proj.mulVec(view.mulVec(vec.extend(vec.Vec3{ 0, 0, 0 }, 1)));
    const ndc_z = clip[2] / clip[3];

    // ── geom: ray/AABB, bezier distance ──
    const V2 = @Vector(2, f32);
    const hit = geom.rayAabb(V2, .{ -5, 0 }, .{ 1, 0 }, .{ -1, -1 }, .{ 1, 1 });
    const p0: V2 = .{ 0, 0 };
    const p1: V2 = .{ 0, 10 };
    const p2: V2 = .{ 10, 10 };
    const p3: V2 = .{ 10, 0 };
    const on_curve = geom.bezierPoint(V2, 0.4, p0, p1, p2, p3);
    const d_on = geom.distPointBezier(V2, 24, on_curve, p0, p1, p2, p3);

    // ── geo: distance/bearing between two cities ──
    const indianapolis = .{ 39.7684, -86.1581 };
    const chicago = .{ 41.8781, -87.6298 };
    const dist_km = geo.distanceKm(f32, indianapolis[0], indianapolis[1], chicago[0], chicago[1]);
    const bearing = geo.bearingDeg(f32, indianapolis[0], indianapolis[1], chicago[0], chicago[1]);
    const compass = geo.cardinal16(f32, bearing);

    // ── grid: 8-neighbor iteration ──
    var cells: [9]u8 = .{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    const g = grid.Grid(u8).initBuffer(&cells, 3, 3);
    var it = g.neighbors8(1, 1);
    var neighbor_sum: usize = 0;
    while (it.next()) |n| neighbor_sum += n.value;

    // ── astro: GMST at a comptime-fixed MJD ──
    const mjd_2026_01_01: f64 = 61041.0; // 2026-01-01 00:00 UTC
    const gmst_rad = try astro.gmst(mjd_2026_01_01);

    try w.print("cross((1,0,0),(0,1,0)) = {d}\n", .{up});
    try w.print("swizzle(1,2,3,4,\"wzyx\") = {d}\n", .{swz});
    try w.print("origin (5 units in front of eye) projects to NDC z = {d:.4}\n", .{ndc_z});
    try w.print("ray from (-5,0) toward AABB[-1,-1]-[1,1]: enters at t = {?d:.2}\n", .{hit});
    try w.print("bezier: point at t=0.4 is {d:.3} from the curve\n", .{d_on});
    try w.print("Indianapolis -> Chicago: {d:.1} km, bearing {d:.1} deg ({s})\n", .{ dist_km, bearing, compass });
    try w.print("sum of 8 neighbors of center cell (value 4) in a 3x3 grid: {d}\n", .{neighbor_sum});
    try w.print("GMST at MJD {d:.1}: {d:.4} rad\n", .{ mjd_2026_01_01, gmst_rad });
    try w.flush();

    const ok = hit != null and d_on < 0.05 and dist_km > 240 and dist_km < 280 and
        std.mem.eql(u8, compass, "NNW") and neighbor_sum == 0 + 1 + 2 + 3 + 5 + 6 + 7 + 8;
    if (!ok) return error.ExampleFailed;
}
