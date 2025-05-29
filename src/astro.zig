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

///   Form a rotation matrix from the Euler angles - three successive
///   rotations about specified Cartesian axes (double precision)
///   Given:
///     ORDER   c*(*)   specifies about which axes the rotations occur
///     PHI     d       1st rotation (radians)
///     THETA   d       2nd rotation (   "   )
///     PSI     d       3rd rotation (   "   )
///   Returned:
///     RMAT    d(3,3)  rotation matrix
///   A rotation is positive when the reference frame rotates
///   anticlockwise as seen looking towards the origin from the
///   positive region of the specified axis.
///   The characters of ORDER define which axes the three successive
///   rotations are about.  A typical value is 'ZXZ', indicating that
///   RMAT is to become the direction cosine matrix corresponding to
///   rotations of the reference frame through PHI radians about the
///   old Z-axis, followed by THETA radians about the resulting X-axis,
///   then PSI radians about the resulting Z-axis.
///   The axis names can be any of the following, in any order or
///   combination:  X, Y, Z, uppercase or lowercase, 1, 2, 3.  Normal
///   axis labelling/numbering conventions apply;  the xyz (=123)
///   triad is right-handed.  Thus, the 'ZXZ' example given above
///   could be written 'zxz' or '313' (or even 'ZxZ' or '3xZ').  ORDER
///   is terminated by length or by the first unrecognized character.
///   Fewer than three rotations are acceptable, in which case the later
///   angle arguments are ignored.  If all rotations are zero, the
///   identity matrix is produced.
pub fn deuler(comptime T: type, order: [3]u8, phi: T, theta: T, psi: T) ![3][3]T {
    var res: [3][3]T = @splat(@splat(1));
    for (order, 0..) |axis, i| {
        var rot: [3][3]T = @splat(@splat(1));
        const angle = switch (i) {
            1 => phi,
            2 => theta,
            else => psi,
        };
        const ang_sin = @sin(angle);
        const ang_cos = @cos(angle);
        switch (axis) {
            'X', 'x', '1' => {
                rot[1][1] = ang_cos;
                rot[1][2] = ang_sin;
                rot[2][1] = -ang_sin;
                rot[2][2] = ang_cos;
            },
            'Y', 'y', '2' => {
                rot[0][0] = ang_cos;
                rot[0][2] = -ang_sin;
                rot[2][0] = ang_sin;
                rot[2][2] = ang_cos;
            },
            'Z', 'z', '3' => {
                rot[0][0] = ang_cos;
                rot[0][1] = ang_sin;
                rot[1][0] = -ang_sin;
                rot[1][1] = ang_cos;
            },
            else => return error.InvalidOrigin,
        }
        res = matmul(T, res, rot);
    }
    return res;
}

fn matmul(comptime T: type, a: [3][3]T, b: [3][3]T) [3][3]T {
    var result: [3][3]f32 = undefined;

    for (0..3) |i| {
        for (0..3) |j| {
            var sum: f32 = 0;
            for (0..3) |k| {
                sum += a[i][k] * b[k][j];
            }
            result[i][j] = sum;
        }
    }

    return result;
}
/// Conversion of Besselian Epoch to Modified Julian Date
/// (double precision)
/// Given:
/// epb      dp       Besselian Epoch
/// The result is the Modified Julian Date (JD - 2400000.5).
/// -
pub fn epb2d(epb: anytype) @TypeOf(epb) {
    return 15019.81352e0 + (epb - 1900e0) * 365.242198781e0;
}

pub fn Gregorian(comptime T: type) type {
    return struct {
        year: T = 0,
        month: T = 0,
        day: T = 0,
        day_fraction: T = 0,
    };
}
/// Modified Julian Date to Gregorian year, month, day,
/// and fraction of a day.
/// Given:
/// DJM      dp     modified Julian Date (JD-2400000.5)
/// Returned:
/// IY       int    year
/// IM       int    month
/// ID       int    day
/// FD       dp     fraction of day
/// J        int    status:
/// 0 = OK
/// -1 = unacceptable date (before 4701BC March 1)
/// The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
/// Last revision:   22 July 2004
pub fn djcl(mjd: anytype) !Gregorian(@TypeOf(mjd)) {
    const T = @TypeOf(mjd);
    var ret: Gregorian(T) = .{};
    if (mjd <= -2395520e0 or mjd >= 1e9) return error.JulianDateOutOfRange;

    ret.day_fraction = @mod(mjd, 1e0);
    if (ret.day_fraction < 0e0)
        ret.day_fraction += 1e0;
    const day_int = math.floor(mjd - ret.day_fraction);
    const jd = day_int + 2400001.0;

    const n4 = 4 * (jd + @divFloor((6 * @divFloor((@divFloor(4 * jd - 17918, 146097)), 4)) + 1, 2) - 37);
    const nd10 = 10 * (@divFloor(@mod(n4 - 237, 1461), 4)) + 5;
    ret.year = @divFloor(n4, 1461) - 4712;
    ret.month = @mod(@divFloor(nd10, 306) + 2, 12) + 1;
    ret.day = @divFloor(@mod(nd10, 306), 10) + 1;
    return ret;
}

pub fn djcal(ndp: anytype, djm: anytype) !Gregorian(@TypeOf(djm)) {
    const T = @TypeOf(djm);
    if ((djm <= -2395520e0) or (djm <= -1e9)) return error.JulianDateOutOfRange;

    // Denominator
    const nfd = math.pow(T, 10, @max(@as(T, ndp), 0));
    const fd = nfd;

    // round date
    const df = djm * fd;

    // separate day and fraction
    var f = @mod(df, fd);
    if (f < 0e0)
        f = f + fd;
    const d = (df - f) / fd;

    // gregorian
    const jd = d + 2400001;
    const n4 = 4 * (jd + @divFloor((6 * @divFloor((@divFloor(4 * jd - 17918, 146097)), 4)) + 1, 2) - 37);
    const nd10 = 10 * (@divFloor(@mod(n4 - 237, 1461), 4)) + 5;

    return .{
        .year = @divFloor(n4, 1461) - 4712,
        .month = @mod(@divFloor(nd10, 306) + 2, 12) + 1,
        .day = @divFloor(@mod(nd10, 306), 10) + 1,
        .day_fraction = f,
    };
}
/// Horizon to equatorial coordinates:  Az,El to HA,Dec
///
/// (double precision)
///
/// Given:
/// AZ      d     azimuth
/// EL      d     elevation
/// PHI     d     observatory latitude
///
/// Returned:
/// HA      d     hour angle
/// DEC     d     declination
///
/// Notes:
///
/// 1)  All the arguments are angles in radians.
///
/// 2)  The sign convention for azimuth is north zero, east +pi/2.
///
/// 3)  HA is returned in the range +/-pi.  Declination is returned
/// in the range +/-pi/2.
///
/// 4)  The latitude is (in principle) geodetic.  In critical
/// applications, corrections for polar motion should be applied.
///
/// 5)  In some applications it will be important to specify the
/// correct type of elevation in order to produce the required
/// type of HA,Dec.  In particular, it may be important to
/// distinguish between the elevation as affected by refraction,
/// which will yield the "observed" HA,Dec, and the elevation
/// in vacuo, which will yield the "topocentric" HA,Dec.  If the
/// effects of diurnal aberration can be neglected, the
/// topocentric HA,Dec may be used as an approximation to the
/// "apparent" HA,Dec.
///
/// 6)  No range checking of arguments is done.
///
/// 7)  In applications which involve many such calculations, rather
/// than calling the present routine it will be more efficient to
/// use inline code, having previously computed fixed terms such
/// as sine and cosine of latitude.
/// -
pub fn de2h(comptime T: type, horizon: Horizon(T)) !Equitorial(T) {
    const sa = @sin(horizon.azimuth);
    const ca = @cos(horizon.azimuth);
    const se = @sin(horizon.elevation);
    const ce = @cos(horizon.elevation);
    const sp = @sin(horizon.obvs_latitude);
    const cp = @cos(horizon.obvs_latitude);

    // HA,Dec as x,y,z
    const x = -ca * ce * sp + se * cp;
    const y = -sa * ce;
    const z = ca * ce * cp + se * sp;

    // To HA,Dec
    const r = @sqrt(x * x + y * y);
    return .{
        .hour_angle = math.atan2(y, x),
        .declanation = math.atan2(z, r),
    };
    // ha = np.where(r == 0, 0, np.arctan2(y, x))
    // dec = np.arctan2(z, r)
    // return ha, dec
}
pub fn Equitorial(comptime T: type) type {
    return struct {
        hour_angle: T,
        declenation: T,
    };
}
pub fn Horizon(comptime T: type) type {
    return struct {
        azimuth: T,
        elevation: T,
        obvs_latitude: T,
    };
}
test {
    std.debug.print("{}\n", .{try clyd(2025, 1, 1)});
    std.debug.print("{}\n", .{try djcl(2460818.372431)});
    std.debug.print("{}\n", .{try djcal(2, @as(f32, 2460818.372431))});
    std.debug.print("{any}\n", .{try deuler(f32, "ZXZ".*, 2, 2, 2)});
}
