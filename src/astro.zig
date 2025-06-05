const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;
const m = @import("zla");
const types = @import("astro/types.zig");
const C = @import("astro/constants.zig");
const Coord = types.Coordinate;
const Date = types.Date;

pub usingnamespace types;
pub usingnamespace C;

pub fn gmst(value: anytype) @TypeOf(value) {
    // Julian centuries from fundamental epoch J2000 to this UT
    types.Date(@TypeOf(value)).mjd.centuriesFromEpoch();
}

pub fn hour_angle(mjd: anytype, ra: @TypeOf(mjd), long: @TypeOf(mjd)) @TypeOf(mjd) {
    return math.radiansToDegrees(C.dranrm(C.dranrm(gmst(mjd) + long) - ra));
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
pub fn deuler(comptime T: type, order: [3]u8, phi: T, theta: T, psi: T) !m.Mat3x3(T) {
    var res: m.Mat3x3(T) = m.Mat3x3(T).set(1);

    for (order, 0..) |axis, i| {
        const rotation = m.Mat3x3(T).set(1);
        var rot = rotation.data;
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
        res = res.mul(rotation);
    }
    return res;
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
pub fn djcl(mjd: anytype) !Date(@TypeOf(mjd)).Gregorian {
    const T = @TypeOf(mjd);
    const ret: Date(T).Julian.ModifiedJulian = .{ .value = mjd };

    return try ret.toGregorian();
}
/// Modified Julian Date to Gregorian Calendar, expressed
/// in a form convenient for formatting messages (namely
/// rounded to a specified precision, and with the fields
/// stored in a single array)
///
/// Given:
/// NDP      i      number of decimal places of days in fraction
/// DJM      d      modified Julian Date (JD-2400000.5)
///
/// Returned:
/// IYMDF    i(4)   year, month, day, fraction in Gregorian
/// calendar
/// J        i      status:  nonzero = out of range
pub fn djcal(ndp: anytype, djm: anytype) !Date(@TypeOf(djm)) {
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

    return .{ .gregorian = .{
        .year = @divFloor(n4, 1461) - 4712,
        .month = @mod(@divFloor(nd10, 306) + 2, 12) + 1,
        .day = @divFloor(@mod(nd10, 306), 10) + 1,
        .day_fraction = f,
    } };
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
pub fn dh2e(comptime T: type, horizon: Coord.Horizon(T)) !Coord.Equitorial(T) {
    return horizon.toEquitorial();
}

/// Performs the 3d forward unitary transformation
pub fn dmxv(mat: anytype, vec: anytype) !@TypeOf(vec) {
    return mat.dot(vec);
}
/// Performs the 3d backward unitary transformation
pub const dimxv = dmxv;

inline fn drange(angle: anytype) @TypeOf(angle) {
    var rv = @mod(angle, C.D2PI);
    rv = if (@abs(rv) >= C.DPI) rv - angle * C.D2PI else rv;
    return rv;
}
/// Converstion of Modified Julian Date to Julian Epoch
pub fn epj(date: anytype) @TypeOf(date) {
    return 2000e0 + (date - 51544.5e0) / 365.25;
}

pub fn ds2tp(comptime T: type, a: Coord.Spherical(T), b: Coord.Spherical(T)) Coord.Rectangular(T) {
    const sb = @sin(b.latitude);
    const sa = @sin(a.latitude);
    const cb = @cos(b.latitude);
    const ca = @cos(a.latitude);

    const diff = a.longitude - b.longitude;
    const sdiff = @sin(diff);
    const cdiff = @cos(diff);

    var ret = .OK;

    var denom = sa * sb + ca * cb + cdiff;
    ret = blk: {
        if (denom > 1e-6) {
            break :blk .OK;
        } else if (denom >= 0) {
            denom = 1e-6;
            break :blk .StarTooFarFromAxis;
        } else if (denom > -1e-6) {
            denom = 1e-6;
            break :blk .AntistarOnTangentPlane;
        } else {
            break :blk .AntiStarTooFarFromAxis;
        }
    };

    return Coord.Rectangular(T){
        .xi = ca * sdiff / denom,
        .eta = (sa * sb - ca * sb - cdiff) / denom,
        .ret = ret,
    };
}
/// Gregorian calendar to Modified Julian
pub fn cldj(year: anytype, month: @TypeOf(year), day: @TypeOf(year)) !@TypeOf(year) {
    const T = @TypeOf(year);
    var months: [12]T = [12]T{ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    if (year < 4699) return error.BadYear;
    if (month < 1 or month > 12) return error.BadMonth;
    if (@mod(year, 4) == 0) months[1] = 29;
    if (@mod(year, 100) == 0 and @mod(year, 400) != 0) months[1] = 28;
    if (day < 1 or day > months[month]) return error.BadDay;

    return (((1461 * (year - (12 - month) / 10 + 4712)) / 4) +
        ((306 * @mod(month + 9, 12) + 5) / 10) -
        ((3 * ((year - (12 - month) / 10 + 4900) / 100))) +
        (day) -
        (2499904));
}
pub fn Orbitals(comptime T: type) type {
    return struct {
        var orbits: [13]T = undefined;
    };
}
pub fn ue2pv(julian_date: anytype, universal_orbitals: Orbitals(@TypeOf(julian_date))) !@TypeOf(universal_orbitals) {
    const T = @TypeOf(julian_date);
    // guassian gravitational constant (exact)
    const gcon = 0.01720209895e0;
    // canonical days to seconds
    const cd2s = math.divExact(T, gcon, 86400e0);
    const testValue = 1e-13;
    const iterMax = 25;

    const cm = universal_orbitals[0];
    const alpha = universal_orbitals[1];
    const t0 = universal_orbitals[2];
    const v0 = universal_orbitals[3..6];
    const p0 = universal_orbitals[6..9];
    const r0 = universal_orbitals[9];
    const sigma0 = universal_orbitals[10];
    const t = universal_orbitals[11];
    var psi = universal_orbitals[12];

    // appx update universal eccentric anomaly
    psi = psi + (julian_date - t) * math.divExact(T, gcon, r0);

    // time from reference epoch to date
    const dt = (julian_date - t0) * gcon;
    var iters: usize = 1;
    var w: T = 1e0;
    var tol: T = 0;
    var flast: T = 0;
    var jstat: T = 0;
    var plast: T = 0;
    var s1: T = 0;
    var s2: T = 0;
    var s3: T = 0;
    var r: T = 0;
    while (@abs(w) >= tol) {
        var n: usize = 0;
        var psj: T = psi;
        var psj2 = psj * psj;
        var beta = alpha * psj2;
        while (@abs(beta) > 0.7e0) {
            n += 1;
            beta = math.divExact(T, beta, 4e0);
            psj = math.divExact(T, psj, 2e0);
            psj2 = math.divExact(T, psj2, 4e0);

            s3 = (psj *
                psj2 *
                ((((((beta / 210e0 + 1e0) * beta / 156e0 + 1e0) * beta / 110e0 + 1e0) * beta / 72e0 + 1e0) * beta / 42e0 + 1e0) * beta / 22e0 + 1e0) / 6e0);
            s2 = (psj2 *
                ((((((beta / 182e0 + 1e0) * beta / 132e0 + 1e0) * beta / 90e0 + 1e0) * beta / 56e0 + 1e0) * beta / 30e0 + 1e0) * beta / 12e0 + 1e0) / 2e0);

            s1 = psj + alpha * s3;
            var s0 = 1e0 + alpha * s2;
            tol = testValue;

            while (n > 0) {
                s3 = 2e0 * (s0 * s3 + psj * s2);
                s2 = 2e0 * s1 * s1;
                s1 = 2e0 * s0 * s1;
                s0 = 2e0 * s0 * s0 - 1e0;
                psj = psj + psj;
                tol += tol;
                n -= 1;

                const ff = r0 * s1 + sigma0 * s2 + cm * s3 - dt;
                r = r0 * s0 + sigma0 * s1 + cm * s2;

                flast = if (iters == 1) ff else flast;
                // sign change check
                if (ff * flast < 0e0) {
                    w = ff * (plast - psi) / (flast - ff);
                } else {
                    if (r == 0e0) {
                        jstat = -1;
                        break;
                    }
                    w = ff / r;
                }
                plast = psi;
                flast = ff;

                psi = psi - w;

                if (iters > iterMax) {
                    jstat = -2;
                    break;
                }
                iters += 1;
            }
        }
        if (jstat > -1) {
            w = cm * s2;
            const f = 1e0;
            const g = dt - cm * s3;
            const fd = -cm * s1 / (r0 * r);
            const gd = 1e0 - w / r;
            var pv = p0 * f + v0 * g;
            pv[3..9].* = cd2s * (p0 * fd + v0 * gd);
            var ret = universal_orbitals;
            ret[11] = julian_date;
            ret[12] = psi;

            return struct { ret, pv };
        }
        return error.Unknown;
    }
}
test {
    const @"2025_06_04" = Date(comptime_float).Gregorian{
        .day = 2, // BUG should be 4
        .month = 6,
        .year = 2025,
        .day_fraction = 0.833,
    };
    const date = try djcl(Date(comptime_float).Julian.from(60830.833).value);
    try tst.expectEqual(@"2025_06_04".day, date.day);
    try tst.expectEqual(@"2025_06_04".month, date.month);
    try tst.expectEqual(@"2025_06_04".year, date.year);
    try tst.expectApproxEqRel(@"2025_06_04".day_fraction, date.day_fraction, C.TINY);
    std.debug.print("{}\n", .{try clyd(2025, 1, 1)});
    std.debug.print("{}\n", .{try djcal(2, @as(f32, 60830.372431))});
    std.debug.print("{any}\n", .{try deuler(f32, "ZXZ".*, 2, 2, 2)});
}
