const std = @import("std");
const Array = std.ArrayList;
const tst = std.testing;
const math = std.math;

/// Pi
pub const DPI = math.pi;

pub const D2PI = math.pi * 2;

/// Seconds to Radians
pub const DS2R = DPI / (12.0 * 3600.0);

/// Acseconds to Radians
pub const AS2R = 4.848136811095359935899141e-6;

/// Epislon
pub const TINY = math.floatEps(f32);

/// Light time for 1 AU (sec)
pub const CR = 499.004782e0;
/// Gravitational radius of the Sun x 2 (2*mu/c**2, AU)
pub const GR2 = 1.974126e-8;
/// B1950
pub const B1950 = 1949.9997904423e0;
/// Degrees to radians
pub const DD2R = 1.745329251994329576923691e-2;
/// Arc seconds in a full circle
pub const TURNAS = 1296000e0;
/// Reference epoch (J2000), MJD
pub const DJM0 = 51544.5e0;
/// Days per Julian century
pub const DJC = 36525e0;
/// Mean sidereal rate (at J2000) in radians per (UT1) second
pub const SR = 7.292115855306589e-5;
/// Earth equatorial radius (metres)
pub const A0 = 6378140e0;
/// Reference spheroid flattening factor and useful function
pub const SPHF = 1e0 / 298.257e0;
pub const SPHB = (1e0 - SPHF) * (1e0 - SPHF);
/// Astronomical unit in metres
pub const AU = 1.49597870e11;

inline fn dmod(A: anytype, B: @TypeOf(A)) !@TypeOf(A) {
    return math.mod(@TypeOf(A), A, B);
}

inline fn dranrm(value: anytype) !@TypeOf(value) {
    return dmod(value, D2PI);
}

pub fn gmst(value: anytype) !@TypeOf(value) {
    // Julian centuries from fundamental epoch J2000 to this UT

    const tu = (value - 51544.5) / 36525.0;
    return try dranrm(try dmod(value, 1) * D2PI + (24110.54841 + (8640184.812866 + (0.093104 - 6.2e-6 * tu) * tu) * tu) + DS2R);
}

pub fn hour_angle(mjd: anytype, ra: @TypeOf(mjd), long: @TypeOf(mjd)) !@TypeOf(mjd) {
    return math.radiansToDegrees(try dranrm(try dranrm(try gmst(mjd) + long) - ra));
}

/// Converts a Gregorian date to (year, day-of-year) in a Julian calendar
/// realigned to match Gregorian dates between 1900-03-01 and 2100-02-28;
/// outside that range the two calendars drift by a day per non-leap century.
pub fn clyd(year: anytype, month: @TypeOf(year), day: @TypeOf(year)) !@Tuple(&[_]type{ @TypeOf(year), @TypeOf(year) }) {
    const T = @TypeOf(year);
    var ret_year: T = 0;
    var ret_day: T = 0;

    if (year >= -4711) {
        if (1 <= month and month <= 12) {
            var month_lengths: [12]T = .{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

            if (@mod(year, 4) == 0 and (@mod(year, 100) != 0 or @mod(year, 400) == 0))
                month_lengths[1] = 29;

            if (day < 1 or day > month_lengths[@as(usize, @intCast(month)) - 1]) return error.BadDay;

            var i = (14 - month) / 12;
            var k = year - i;
            var j = (1461 * (k + 4800) / 4 + (367 * (month - 2 + 12 * i)) / 12 - (3 * ((k + 4900) / 100)) / 4 + day - 3660);
            k = (j - 1) / 1461;
            const l = j - 1461 * k;
            const n = (l - 1) / 365 - l / 1461;
            j = ((80 * (l - 365 * n + 30)) / 2447) / 11;
            i = n + j;

            ret_day = 59 + l - 365 * i + ((4 - n) / 4) * (1 - j);
            ret_year = 4 * k + i - 4716;
        } else return error.BadMonth;
    } else return error.BadYear;

    return .{
        ret_year,
        ret_day,
    };
}
pub fn Cartesian(comptime T: type) type {
    return struct {
        x: T,
        y: T,
        z: T,
    };
}
pub fn Spherical(comptime T: type) type {
    return struct {
        latitude: T,
        longitude: T,
    };
}

/// Longitude is +ve anticlockwise looking from the +ve latitude pole; the
/// x axis is at zero longitude/latitude, z axis at the +ve latitude pole.
/// At either pole, longitude is returned as zero.
pub fn dcc2s(comptime T: type, coord: Cartesian(T)) Spherical(T) {
    const r = @sqrt(coord.x * coord.x + coord.y * coord.y);
    return .{
        .latitude = math.atan2(coord.z, r),
        .longitude = math.atan2(coord.y, coord.x),
    };
}

/// Inverse of `dcc2s`: same longitude/latitude convention (see there).
pub fn dcs2c(comptime T: type, sphere: Spherical(T)) Cartesian(T) {
    const right_acention = sphere.longitude;
    const declanation = sphere.latitude;
    return .{
        .x = @cos(right_acention) * @cos(declanation),
        .y = @sin(right_acention) * @cos(declanation),
        .z = @sin(declanation),
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

test gmst {
    // MJD 51544.5 (J2000 epoch) should give gmst in radians
    const result = try gmst(51544.5);
    try tst.expect(result > 4.0 and result < 6.0);
}

test "dcs2c/dcc2s round trip" {
    const c1 = dcs2c(f64, .{ .longitude = 0.5, .latitude = 0.3 });
    const s2 = dcc2s(f64, c1);
    try tst.expectApproxEqAbs(@as(f64, 0.5), s2.longitude, 1e-10);
    try tst.expectApproxEqAbs(@as(f64, 0.3), s2.latitude, 1e-10);
}

test "clyd known date" {
    const result = try clyd(2025, 1, 1);
    // clyd converts to Julian calendar; 2025/1/1 Gregorian -> 2098/338 Julian
    try tst.expect(result[0] == 2098);
    try tst.expect(result[1] == 338);
}
