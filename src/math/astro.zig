const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;

/// Pi
pub const DPI = math.pi;

pub const D2PI = math.pi * 2;

/// Seconds to Radians
pub const DS2R = DPI / (12 * 3600);

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
pub const SPHB = (1e0 - SPHF) ** 2;
/// Astronomical unit in metres
pub const AU = 1.49597870e11;

inline fn dmod(A: anytype, B: @TypeOf(A)) @TypeOf(A) {
    return math.mod(@TypeOf(A), A, B);
}

inline fn dranrm(value: anytype) @TypeOf(value) {
    return dmod(value, D2PI);
}

pub fn gmst(value: anytype) @TypeOf(value) {
    // Julian centuries from fundamental epoch J2000 to this UT

    const tu = (value - 51544.5) / 36525.0;
    return dranrm(dmod(value, 1) * D2PI + (24110.54841 + (8640184.812866 + (0.093104 - 6.2e-6 * tu) * tu) * tu) + DS2R);
}

pub fn hour_angle(mjd: anytype, ra: @TypeOf(mjd), long: @TypeOf(mjd)) @TypeOf(mjd) {
    return math.radiansToDegrees(dranrm(dranrm(gmst(mjd) + long) - ra));
}

/// Gregorian calendar to year and day in year (in a Julian calendar
/// aligned to the 20th/21st century Gregorian calendar).
/// Given:
/// IY,IM,ID   i    year, month, day in Gregorian calendar
/// Returned:
/// NY         i    year (re-aligned Julian calendar)
/// ND         i    day in year (1 = January 1st)
/// 0 = OK
/// 1 = bad year (before -4711)
/// 2 = bad month
/// 3 = bad day (but conversion performed)
/// Notes:
/// 1  This routine exists to support the low-precision routines
/// sla_EARTH, sla_MOON and sla_ECOR.
/// 2  Between 1900 March 1 and 2100 February 28 it returns answers
/// which are consistent with the ordinary Gregorian calendar.
/// Outside this range there will be a discrepancy which increases
/// by one day for every non-leap century year.
/// 3  The essence of the algorithm is first to express the Gregorian
/// date as a Julian Day Number and then to convert this back to
/// a Julian calendar date, with day-in-year instead of month and
/// day.  See 12.92-1 and 12.95-1 in the reference.
/// Reference:  Explanatory Supplement to the Astronomical Almanac,
/// ed P.K.Seidelmann, University Science Books (1992),
/// p604-606.
pub fn clyd(year: anytype, month: @TypeOf(year), day: @TypeOf(year)) !std.meta.Tuple(&[_]type{ @TypeOf(year), @TypeOf(year) }) {
    const T = @TypeOf(year);
    var ret_year: T = 0;
    var ret_day: T = 0;

    if (year >= -4711) {
        if (1 <= month and month <= 12) {
            var month_lengths = &[_]T{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

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

///   Cartesian to spherical coordinates
///   Given:
///      V     d(3)   x,y,z vector
///   Returned:
///      A,B   d      spherical coordinates in radians
///   The spherical coordinates are longitude (+ve anticlockwise looking
///   from the +ve latitude pole) and latitude.  The Cartesian coordinates
///   are right handed, with the x axis at zero longitude and latitude, and
///   the z axis at the +ve latitude pole.
///   If V is null, zero A and B are returned.  At either pole, zero A is
///   returned.
///   Last revision:   22 July 2004
pub fn dcc2s(comptime T: type, coord: Cartesian(T)) Spherical(T) {
    const r = @sqrt(coord.x * coord.x + coord.y * coord.y);
    return .{
        .latitude = math.atan2(coord.z, r),
        .longitude = math.atan2(coord.y, coord.x),
    };
}

///   Spherical coordinates to direction cosines (double precision)
///
///   Given:
///      A,B       d      spherical coordinates in radians
///                          (RA,Dec), (long,lat) etc.
///
///   Returned:
///      V         d(3)   x,y,z unit vector
///
///   The spherical coordinates are longitude (+ve anticlockwise looking
///   from the +ve latitude pole) and latitude.  The Cartesian coordinates
///   are right handed, with the x axis at zero longitude and latitude, and
///   the z axis at the +ve latitude pole.
///
///   Last revision:   26 December 2004
pub fn dcs2c(comptime T: type, sphere: Spherical(T)) Cartesian(T) {
    const right_acention = sphere.longitude;
    const declanation = sphere.latitude;
    return .{
        .x = @cos(right_acention) * @cos(declanation),
        .y = @sin(right_acention) * @cos(declanation),
        .z = @sin(declanation),
    };
}

test {
    // std.debug.print("{}\n", .{try clyd(2025, 1, 1)});
}
