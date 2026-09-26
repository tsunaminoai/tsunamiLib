const std = @import("std");

/// Mean earth radius [km].
pub fn earthRadiusKm(comptime T: type) T {
    return 6371.0;
}

/// Effective earth radius [km] for the standard 4/3-earth refraction model
/// used by weather radar beam propagation (Doviak & Zrnić 1993, §2.2.3).
pub fn fourThirdsRadiusKm(comptime T: type) T {
    return 4.0 / 3.0 * earthRadiusKm(T);
}

/// Great-circle distance [km] between two lat/lon points [deg]. Haversine,
/// computed in f64 internally regardless of `T` for numerical headroom, then
/// cast back.
pub fn distanceKm(comptime T: type, lat_a: T, lon_a: T, lat_b: T, lon_b: T) T {
    const rad = std.math.pi / 180.0;
    const la1: f64 = @as(f64, lat_a) * rad;
    const la2: f64 = @as(f64, lat_b) * rad;
    const dla: f64 = la2 - la1;
    const dlo: f64 = (@as(f64, lon_b) - @as(f64, lon_a)) * rad;
    const sdla = @sin(dla / 2.0);
    const sdlo = @sin(dlo / 2.0);
    const h = sdla * sdla + @cos(la1) * @cos(la2) * sdlo * sdlo;
    const cl = std.math.clamp(h, 0.0, 1.0);
    return @floatCast(2.0 * earthRadiusKm(f64) * std.math.asin(@sqrt(cl)));
}

/// Initial compass bearing [deg, 0=north, 90=east] from (lat_a,lon_a) to (lat_b,lon_b).
pub fn bearingDeg(comptime T: type, lat_a: T, lon_a: T, lat_b: T, lon_b: T) T {
    const rad = std.math.pi / 180.0;
    const la1 = lat_a * rad;
    const la2 = lat_b * rad;
    const dlo = (lon_b - lon_a) * rad;
    const y = @sin(dlo) * @cos(la2);
    const x = @cos(la1) * @sin(la2) - @sin(la1) * @cos(la2) * @cos(dlo);
    const b = std.math.radiansToDegrees(std.math.atan2(y, x));
    return if (b < 0) b + 360.0 else b;
}

/// Destination point [deg] `dist_km` along initial `bearing_deg` from (lat,lon).
/// Returns `.{ lat, lon }`.
pub fn destination(comptime T: type, lat: T, lon: T, bearing_deg: T, dist_km: T) [2]T {
    const rad = std.math.pi / 180.0;
    const ang = dist_km / earthRadiusKm(T);
    const la1 = lat * rad;
    const lo1 = lon * rad;
    const brg = bearing_deg * rad;
    const la2 = std.math.asin(@sin(la1) * @cos(ang) + @cos(la1) * @sin(ang) * @cos(brg));
    const lo2 = lo1 + std.math.atan2(
        @sin(brg) * @sin(ang) * @cos(la1),
        @cos(ang) - @sin(la1) * @sin(la2),
    );
    return .{ std.math.radiansToDegrees(la2), std.math.radiansToDegrees(lo2) };
}

/// Project a lat/lon to km east / km north of `ref_lat`/`ref_lon`.
/// Equirectangular with a cos(lat) longitude factor, evaluated at the
/// reference latitude. Good to well under a km across a domain of a few
/// hundred km. Returns `.{ east_km, north_km }`.
pub fn enuKm(comptime T: type, ref_lat: T, ref_lon: T, lat: T, lon: T) [2]T {
    const rad: T = std.math.pi / 180.0;
    const east = (lon - ref_lon) * kmPerDegLon(T) * @cos(ref_lat * rad);
    const north = (lat - ref_lat) * kmPerDegLat(T);
    return .{ east, north };
}

/// Inverse of `enuKm`: recover lat/lon from km east/north of `ref_lat`/`ref_lon`.
pub fn enuKmInverse(comptime T: type, ref_lat: T, ref_lon: T, east_km: T, north_km: T) [2]T {
    const rad: T = std.math.pi / 180.0;
    const lat = ref_lat + north_km / kmPerDegLat(T);
    const lon = ref_lon + east_km / (kmPerDegLon(T) * @cos(ref_lat * rad));
    return .{ lat, lon };
}

fn kmPerDegLat(comptime T: type) T {
    return 110.574;
}
fn kmPerDegLon(comptime T: type) T {
    return 111.320;
}

const CARDINALS = [16][:0]const u8{
    "N", "NNE", "NE", "ENE",
    "E", "ESE", "SE", "SSE",
    "S", "SSW", "SW", "WSW",
    "W", "WNW", "NW", "NNW",
};

/// Nearest 16-point compass label for a bearing [deg].
pub fn cardinal16(comptime T: type, bearing_deg: T) [:0]const u8 {
    const sector: usize = @intFromFloat(@round(@mod(bearing_deg, 360.0) / 22.5));
    return CARDINALS[sector % 16];
}

/// Where a radar gate physically sits under the 4/3-earth model: height above
/// the antenna and great-circle ground range from it, given slant range and
/// launch elevation.
pub const BeamPoint = struct { height_km: f32, ground_km: f32 };

/// The ray that reaches a point: slant range and launch elevation.
pub const BeamRay = struct { range_km: f32, el_rad: f32 };

/// Height/ground-range of a gate at (`range_km`, `el_rad`) under earth radius
/// `ke_radius_km` (typically `fourThirdsRadiusKm`). Numerically stable form:
/// see module docs in the port source — avoids catastrophic cancellation of
/// the ~7e7 km² earth-radius term.
pub fn beamForward(ke_radius_km: f32, range_km: f32, el_rad: f32) BeamPoint {
    const r = range_km;
    const big_r = ke_radius_km;
    const sin_el = @sin(el_rad);
    const x = r * @cos(el_rad);
    const z = big_r + r * sin_el;
    const a = r * r + 2.0 * r * big_r * sin_el;
    const h = a / (@sqrt(big_r * big_r + a) + big_r);
    return .{ .height_km = h, .ground_km = big_r * std.math.atan2(x, z) };
}

/// Exact inverse of `beamForward`: which ray reaches (`ground_km`, `height_km`).
pub fn beamInverse(ke_radius_km: f32, ground_km: f32, height_km: f32) BeamRay {
    const big_r = ke_radius_km;
    const h = height_km;
    const theta = ground_km / big_r;
    const rho = big_r + h;
    const sh = @sin(0.5 * theta);
    const sh2 = sh * sh;
    const r = @sqrt(h * h + 4.0 * big_r * rho * sh2);
    const el = std.math.atan2(h * @cos(theta) - 2.0 * big_r * sh2, rho * @sin(theta));
    return .{ .range_km = r, .el_rad = el };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "distanceKm: zero for identical points, sane for a known pair" {
    try testing.expectApproxEqAbs(@as(f32, 0), distanceKm(f32, 39.7, -86.28, 39.7, -86.28), 1e-4);
    // Indianapolis to Chicago, roughly 260 km.
    const d = distanceKm(f32, 39.708, -86.28, 41.604, -88.085);
    try testing.expect(d > 240 and d < 280);
}

test "bearingDeg cardinal directions" {
    try testing.expectApproxEqAbs(@as(f32, 0.0), bearingDeg(f32, 0, 0, 1, 0), 1e-3);
    try testing.expectApproxEqAbs(@as(f32, 90.0), bearingDeg(f32, 0, 0, 0, 1), 1e-3);
    try testing.expectApproxEqAbs(@as(f32, 180.0), bearingDeg(f32, 0, 0, -1, 0), 1e-3);
    try testing.expectApproxEqAbs(@as(f32, 270.0), bearingDeg(f32, 0, 0, 0, -1), 1e-3);
}

test "destination inverts distanceKm/bearingDeg round trip" {
    const lat: f64 = 39.7;
    const lon: f64 = -86.28;
    const brg: f64 = 47.0;
    const dist: f64 = 120.0;
    const dst = destination(f64, lat, lon, brg, dist);
    const back_d = distanceKm(f64, lat, lon, dst[0], dst[1]);
    const back_b = bearingDeg(f64, lat, lon, dst[0], dst[1]);
    try testing.expectApproxEqAbs(dist, back_d, 1e-2);
    try testing.expectApproxEqAbs(brg, back_b, 1e-2);
}

test "enuKm and its inverse round-trip" {
    const ref_lat: f32 = 39.708;
    const ref_lon: f32 = -86.28;
    const lat: f32 = 40.1;
    const lon: f32 = -85.9;
    const enu = enuKm(f32, ref_lat, ref_lon, lat, lon);
    const back = enuKmInverse(f32, ref_lat, ref_lon, enu[0], enu[1]);
    try testing.expectApproxEqAbs(lat, back[0], 1e-4);
    try testing.expectApproxEqAbs(lon, back[1], 1e-4);
}

test "cardinal16 bands" {
    try testing.expectEqualStrings("N", cardinal16(f32, 0));
    try testing.expectEqualStrings("N", cardinal16(f32, 11.0));
    try testing.expectEqualStrings("NNE", cardinal16(f32, 11.3));
    try testing.expectEqualStrings("E", cardinal16(f32, 90));
    try testing.expectEqualStrings("NNW", cardinal16(f32, 337.5));
    try testing.expectEqualStrings("N", cardinal16(f32, 354.0));
}

test "beam forward/inverse round trip across the scan domain" {
    const ke = fourThirdsRadiusKm(f32);
    const els = [_]f32{ -0.5, 0.0, 0.5, 1.3, 2.4, 4.0, 6.4, 10.0, 19.5, 45.0 };
    var r: f32 = 0.25;
    while (r <= 460.0) : (r *= 1.7) {
        for (els) |el_deg| {
            const el = std.math.degreesToRadians(el_deg);
            const p = beamForward(ke, r, el);
            const back = beamInverse(ke, p.ground_km, p.height_km);
            try testing.expectApproxEqAbs(r, back.range_km, 1e-4);
            try testing.expectApproxEqAbs(el, back.el_rad, 1e-4);
        }
    }
}

test "beam inverse survives near-radar f32 cancellation" {
    const ke = fourThirdsRadiusKm(f32);
    const cases = [_]struct { s: f32, h: f32, r: f32 }{
        .{ .s = 1.0, .h = 0.5, .r = 1.11806 },
        .{ .s = 2.0, .h = 0.5, .r = 2.06161 },
        .{ .s = 150.0, .h = 5.0, .r = 150.12548 },
    };
    for (cases) |c| {
        try testing.expectApproxEqAbs(c.r, beamInverse(ke, c.s, c.h).range_km, 1e-4);
    }
}

test "beam: infinite earth radius recovers the flat-earth limit" {
    const flat: f32 = 1.0e9;
    const r: f32 = 100.0;
    const el = std.math.degreesToRadians(@as(f32, 5.0));
    const p = beamForward(flat, r, el);
    try testing.expectApproxEqRel(r * @sin(el) + 1.0, p.height_km + 1.0, 1e-5);
    try testing.expectApproxEqRel(r * @cos(el), p.ground_km, 1e-5);
}
