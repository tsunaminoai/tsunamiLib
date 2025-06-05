const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;
const m = @import("zla");
const C = @import("constants.zig");

pub const DirectionCosines = m.Vec3;

pub const Coordinate = struct {
    pub fn Spherical(comptime T: type) type {
        return struct {
            latitude: T,
            longitude: T,

            pub usingnamespace m.Vec2;
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
            pub fn direction_cosines(self: @This()) Cartesian(T) {
                return .{
                    .x = @cos(self.right_acention) * @cos(self.declanation),
                    .y = @sin(self.right_acention) * @cos(self.declanation),
                    .z = @sin(self.declanation),
                };
            }
        };
    }

    pub fn Cartesian(comptime T: type) type {
        return struct {
            // x: T = 0,
            // y: T = 0,
            // z: T = 0,

            pub usingnamespace m.Vec3;

            const Self = @This();
            /// Position angle of one celestial direction with respect to another.
            /// # The result is the bearing (position angle), in radians, of point
            /// # V2 with respect to point V1.  It is in the range +/- pi.  The
            /// # sense is such that if V2 is a small distance east of V1, the
            /// # bearing is about +pi/2.  Zero is returned if the two points
            /// # are coincident.
            pub fn bearingAngle(self: Self, other: Self) T {
                var v1: Self = self;
                const w = @sqrt(self.x * self.x + self.y * self.y + self.z * self.z);
                if (w != 0e0) {
                    v1.x /= w;
                    v1.y /= w;
                    v1.z /= w;
                }
                const sq = other.y * v1.x - other.x * v1.y;
                const cq = other.z * (v1.x * v1.x + v1.y * v1.y) - v1.z * (other.x * v1.x + other.y * v1.y);
                return std.math.atan2(sq, if (cq == 0 and sq == 0) 1 else cq);
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
            pub fn toSpherical(self: Self) Spherical(T) {
                const r = self.lengthSq();
                return .{
                    .latitude = math.atan2(self.x(), r),
                    .longitude = math.atan2(self.y(), self.x()),
                };
            }
        };
    }

    pub fn Rectangular(comptime T: type) type {
        return struct {
            ret: Ret,
            xi: T,
            eta: T,

            pub const Ret = enum {
                OK,
                StarTooFarFromAxis,
                AntistarOnTangentPlane,
                AntistarTooFarFromAxis,
            };
        };
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
            pub fn toEquitorial(self: @This()) Equitorial(T) {
                const sa = @sin(self.azimuth);
                const ca = @cos(self.azimuth);
                const se = @sin(self.elevation);
                const ce = @cos(self.elevation);
                const sl = @sin(self.obsv_latitude);
                const cl = @cos(self.obvs_latitude);

                const x = -ca * ce * sl + se * cl;
                const y = -sa * ce;
                const z = ca * ce * cl + se * sl;

                const r = @sqrt(x * x + y * y);
                return .{
                    .hour_angle = math.atan2(y, x),
                    .declenation = math.atan2(z, r),
                };
            }
            const SineCosine = struct { sa: T, ca: T, se: T, ce: T, sp: T, cp: T };
            pub fn sine_cosine(self: @This()) SineCosine {
                return .{
                    .sa = @sin(self.azimuth),
                    .ca = @cos(self.azimuth),
                    .se = @sin(self.elevation),
                    .ce = @cos(self.elevation),
                    .sp = @sin(self.obvs_latitude),
                    .cp = @cos(self.obvs_latitude),
                };
            }
        };
    }
};

pub fn Date(comptime T: type) type {
    return union(enum) {
        gregorian: Gregorian,
        julian: Julian,
        mjd: Julian.ModifiedJulian,

        pub const Gregorian = struct {
            year: T = 0,
            month: T = 0,
            day: T = 0,
            day_fraction: T = 0,
        };

        pub const Julian = struct {
            value: T = 0,
            pub fn from(val: anytype) Julian {
                return Date(@TypeOf(val)).Julian{ .value = val };
            }
            pub fn toModifiedJulian(self: Julian) ModifiedJulian {
                return .{ .value = self.value - 2400000.5 };
            }
            pub const ModifiedJulian = struct {
                value: T = 0,
                pub fn centuriesFromEpoch(self: ModifiedJulian) ModifiedJulian {
                    const tu = (self.value - 51544.5) / 36525.0;
                    return .{ .value = C.dranrm(C.dmod(self.value, 1) * C.D2PI + (24110.54841 + (8640184.812866 + (0.093104 - 6.2e-6 * tu) * tu) * tu) + C.DS2R) };
                }
                pub fn toJulian(self: ModifiedJulian) Julian {
                    return .{ .value = self.value + 2400000.5 };
                }
                pub fn from(val: anytype) ModifiedJulian {
                    return Date(@TypeOf(val)).Julian.ModifiedJulian{ .value = val };
                }

                pub fn toGregorian(self: ModifiedJulian) !Gregorian {
                    var ret: Date(T).Gregorian = .{};
                    if (self.value <= -2395520e0 or self.value >= 1e9) return error.DateOutOfRange;

                    ret.day_fraction = self.value - @floor(self.value);
                    if (ret.day_fraction < 0e0)
                        ret.day_fraction += 1e0;
                    const day_int = math.floor(self.value - ret.day_fraction);
                    const jd = ModifiedJulian.from(day_int).toJulian().value;

                    const n4 = 4 * (jd + @divFloor((6 * @divFloor((@divFloor(4 * jd - 17918, 146097)), 4)) + 1, 2) - 37);
                    const nd10 = 10 * (@divFloor(@mod(n4 - 237, 1461), 4)) + 5;
                    ret.year = @divFloor(n4, 1461) - 4712;
                    ret.month = @mod(@divFloor(nd10, 306) + 2, 12) + 1;
                    ret.day = @divFloor(@mod(nd10, 306), 10) + 1;

                    return ret;
                }
            };
        };
    };
}
