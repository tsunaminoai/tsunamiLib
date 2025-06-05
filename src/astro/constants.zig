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

pub inline fn dmod(A: anytype, B: @TypeOf(A)) @TypeOf(A) {
    return math.mod(@TypeOf(A), A, B);
}

pub inline fn dranrm(value: anytype) @TypeOf(value) {
    return dmod(value, D2PI);
}
pub inline fn Cast(comptime To: type, value: anytype) !@TypeOf(To) {
    return switch (@typeInfo(To)) {
        .float, .comptime_float => switch (@typeInfo(@TypeOf(value))) {
            .int, .comptime_int => @intFromFloat(value),
            .float, .comptime_float => @floatCast(value),
            else => error.UnhandledTypeCast,
        },
        .int, .comptime_int => switch (@typeInfo(@TypeOf(value))) {
            .int, .comptime_int => @intCast(value),
            .float, .comptime_float => @floatFromInt(value),
            else => error.UnhandledTypeCast,
        },
        else => error.UnhandledTypeCast,
    };
}
