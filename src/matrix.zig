const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;
pub const Mat = struct { T: type, size_n: usize, ptr: *anyopaque };

pub fn Matrix(comptime T: type, comptime N: usize) type {
    return struct {
        data: [N][N]T,

        const Self = @This();

        /// Initialize matrix with given values
        pub fn init(values: [N][N]T) Self {
            return .{ .data = values };
        }
        pub fn mat(self: *Self) Mat {
            return .{
                .T = T,
                .size_n = N,
                .ptr = self,
            };
        }

        const MatResult = struct {
            inverse: Self,
            determinant: T,
            soltuions: [N]T,
        };
        pub fn dmat(self: *Self, unknowns: *[N]T) !MatResult {
            //TODO: copy self into mat result to work on
            //TODO: Return errors?
            var iw: [N]T = [_]T{0} ** N;
            const sfa = 1e-20;

            var jf: i32 = 0;
            var d: T = 1e0;
            for (0..N) |k| {
                var amx = @abs(self.data[k][k]);
                var imx = k;
                if (k != N) {
                    for (k..N) |i| {
                        const t = @abs(self.data[i][k]);
                        if (t > amx) {
                            amx = t;
                            imx = i;
                        }
                    }
                }
                if (amx < sfa)
                    jf = -1
                else {
                    if (imx != k) {
                        for (0..N) |j| {
                            const t = self.data[k][j];
                            self.data[k][j] = self.data[imx][j];
                            self.data[imx][j] = t;
                        }
                        const t = unknowns[k];
                        unknowns[k] = unknowns[imx];
                        unknowns[imx] = t;
                        d = -d;
                    }

                    iw[k] = imx;
                    var akk = self.data[k][k];
                    d = d * akk;
                    if (@abs(d) < sfa)
                        jf = -1
                    else {
                        akk = 1e0 / akk;
                        self.data[k][k] = akk;
                        for (0..N) |j|
                            self.data[k][j] = if (j != k) self.data[k][j] * akk else self.data[k][j];

                        const yk = unknowns[k] * akk;
                        unknowns[k] = yk;
                        for (0..N) |i| {
                            const aik = self.data[i][k];
                            if (i != k) {
                                for (0..N) |j| {
                                    self.data[i][j] = if (j != k) self.data[i][j] - aik * self.data[k][j] else self.data[i][j];
                                }
                                unknowns[i] = unknowns[i] - aik * yk;
                            }
                        }

                        for (0..N) |i| {
                            self.data[i][k] = if (i != k) -self.data[i][k] * akk else self.data[i][k];
                        }
                    }
                }
            }
            if (jf != 0) {
                d = 0;
            } else {
                for (0..N) |k| {
                    const np1mk = N + 1 - k;
                    const ki = iw[np1mk];
                    if (np1mk != ki) {
                        for (0..N) |i| {
                            const t = self.data[i][np1mk];
                            self.data[i][np1mk] = self.data[i][ki];
                            self.data[i][ki] = t;
                        }
                    }
                }
            }
            return MatResult{
                .determinant = d,
                .inverse = self,
                .soltuions = iw,
            };
        }

        /// Matrix addition (element-wise)
        pub fn add(self: Self, other: Self) Self {
            var result: Self = undefined;
            for (0..N) |i| {
                for (0..N) |j| {
                    result.data[i][j] = self.data[i][j] + other.data[i][j];
                }
            }
            return result;
        }

        /// Matrix subtraction (element-wise)
        pub fn sub(self: Self, other: Self) Self {
            var result: Self = undefined;
            for (0..N) |i| {
                for (0..N) |j| {
                    result.data[i][j] = self.data[i][j] - other.data[i][j];
                }
            }
            return result;
        }

        /// Matrix multiplication
        pub fn multiply(self: Self, other: Self) Self {
            var result: Self = undefined;
            for (0..N) |i| {
                for (0..N) |j| {
                    var sum: T = 0;
                    for (0..N) |k| {
                        sum += self.data[i][k] * other.data[k][j];
                    }
                    result.data[i][j] = sum;
                }
            }
            return result;
        }

        /// Matrix transposition
        pub fn transpose(self: Self) Self {
            var result: Self = undefined;
            for (0..N) |i| {
                for (0..N) |j| {
                    result.data[i][j] = self.data[j][i];
                }
            }
            return result;
        }

        /// Optimized dot product using SIMD and matrix transposition
        pub fn dot(self: @This(), other: @This()) @This() {
            const other_t = other.transpose();
            var result: @This() = undefined;

            inline for (0..N) |i| {
                const row_a = self.data[i];
                inline for (0..N) |j| {
                    const row_b = other_t.data[j];
                    result.data[i][j] = simdDot(row_a, row_b);
                }
            }
            return result;
        }

        fn simdDot(a: [N]T, b: [N]T) T {
            const vec_size = comptime blk: {
                if (@typeInfo(T) == .float) break :blk 4;
                if (N % 4 == 0) break :blk 4;
                if (N % 2 == 0) break :blk 2;
                break :blk 1;
            };

            const num_vecs = N / vec_size;
            var sum: T = 0;

            inline for (0..num_vecs) |k| {
                const offset = k * vec_size;
                const vec_a: @Vector(vec_size, T) = a[offset .. offset + vec_size].*;
                const vec_b: @Vector(vec_size, T) = b[offset .. offset + vec_size].*;
                sum += @reduce(.Add, vec_a * vec_b);
            }

            // Handle remaining elements
            inline for (num_vecs * vec_size..N) |k| {
                sum += a[k] * b[k];
            }
            return sum;
        }

        /// Formatting for printing
        pub fn format(
            self: Self,
            comptime fmt: []const u8,
            options: std.fmt.FormatOptions,
            writer: anytype,
        ) !void {
            _ = fmt;
            _ = options;
            for (self.data) |row| {
                try writer.print("| ", .{});
                for (row) |val| {
                    try writer.print("{d:8.2} ", .{val});
                }
                try writer.print("|\n", .{});
            }
        }
    };
}

// Example usage
test "Matrix operations" {
    const M = Matrix(f32, 2);

    const a = M.init(.{
        .{ 1, 2 },
        .{ 3, 4 },
    });

    const b = M.init(.{
        .{ 5, 6 },
        .{ 7, 8 },
    });

    const sum = a.add(b);
    const product = a.multiply(b);
    const transposed = a.transpose();

    try std.testing.expectEqual(sum.data[0][0], 6);
    try std.testing.expectEqual(product.data[1][1], 50);
    try std.testing.expectEqual(transposed.data[0][1], 3);
}

// Assume the Matrix struct is defined as in previous responses.

test "Matrix basic operations and dot product" {
    const M2 = Matrix(f32, 2);
    const a2 = M2.init(.{
        .{ 1, 2 },
        .{ 3, 4 },
    });
    const b2 = M2.init(.{
        .{ 5, 6 },
        .{ 7, 8 },
    });

    // Addition
    const sum2 = a2.add(b2);
    try std.testing.expectEqual(sum2.data[0][0], 6);
    try std.testing.expectEqual(sum2.data[0][1], 8);
    try std.testing.expectEqual(sum2.data[1][0], 10);
    try std.testing.expectEqual(sum2.data[1][1], 12);

    // Subtraction
    const diff2 = a2.sub(b2);
    try std.testing.expectEqual(diff2.data[0][0], -4);
    try std.testing.expectEqual(diff2.data[0][1], -4);
    try std.testing.expectEqual(diff2.data[1][0], -4);
    try std.testing.expectEqual(diff2.data[1][1], -4);

    // Multiplication
    const prod2 = a2.multiply(b2);
    try std.testing.expectEqual(prod2.data[0][0], 19); // 1*5 + 2*7
    try std.testing.expectEqual(prod2.data[0][1], 22); // 1*6 + 2*8
    try std.testing.expectEqual(prod2.data[1][0], 43); // 3*5 + 4*7
    try std.testing.expectEqual(prod2.data[1][1], 50); // 3*6 + 4*8

    // Transpose
    const trans2 = a2.transpose();
    try std.testing.expectEqual(trans2.data[0][1], 3);
    try std.testing.expectEqual(trans2.data[1][0], 2);

    // Dot product (should be the same as multiply for matrices)
    const dot2 = a2.dot(b2);
    try std.testing.expectEqual(dot2.data[0][0], 19);
    try std.testing.expectEqual(dot2.data[0][1], 22);
    try std.testing.expectEqual(dot2.data[1][0], 43);
    try std.testing.expectEqual(dot2.data[1][1], 50);

    // 4x4 test for dot product
    const M4 = Matrix(f32, 4);
    const a4 = M4.init(.{
        .{ 1, 2, 3, 4 },
        .{ 5, 6, 7, 8 },
        .{ 9, 10, 11, 12 },
        .{ 13, 14, 15, 16 },
    });
    const b4 = M4.init(.{
        .{ 16, 15, 14, 13 },
        .{ 12, 11, 10, 9 },
        .{ 8, 7, 6, 5 },
        .{ 4, 3, 2, 1 },
    });
    const dot4 = a4.dot(b4);
    // Spot-check a few values
    try std.testing.expectEqual(dot4.data[0][0], 80.0); // 1*16 + 2*12 + 3*8 + 4*4
    try std.testing.expectEqual(dot4.data[3][3], 386.0); // 13*13 + 14*9 + 15*5 + 16*1
    try std.testing.expectEqual(dot4.data[1][2], 188.0); // 5*14 + 6*10 + 7*6 + 8*2
}
