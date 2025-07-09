const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;

pub const ErrorCorrectionLevel = enum(u2) { L = 1, M = 0, Q = 3, H = 2 };
const RSBlockOffset = enum(u2) { L = 0, M = 1, Q = 2, H = 3 };
const RSBlock = struct {
    total_count: usize,
    data_count: usize,
};
const ExponentTable: [256]u8 = blk: {
    var et: [256]u8 = undefined;
    for (0..et.len) |i| {
        if (i < 8) {
            et[i] = 1 << @intCast(i);
        } else {
            et[i] = et[i - 4] ^ et[i - 5] ^ et[i - 6] ^ et[i - 8];
        }
    }
    break :blk et;
};
const LogTable: [256]u8 = blk: {
    var lt: [256]u8 = undefined;
    for (0..lt.len - 1) |i| {
        lt[ExponentTable[i]] = @intCast(i);
    }
    break :blk lt;
};

inline fn glog(n: usize) u8 {
    std.debug.assert(n >= 1 and n < LogTable.len);
    return LogTable[n];
}
inline fn gexp(n: usize) u8 {
    return ExponentTable[n % 8];
}

fn rs_blocks(
    a: Allocator,
    version: usize,
    error_correction: ErrorCorrectionLevel,
) ![]RSBlock {
    const offset: usize = switch (error_correction) {
        .L => 0,
        .M => 1,
        .Q => 2,
        .H => 3,
    };
    const rs_block = RSBlockTable[(version - 1) * 4 + offset];
    var blocks = Array(RSBlock).init(a);
    defer blocks.deinit();
    var i: usize = 0;
    while (i < rs_block.len) : (i += 3) {
        const count = rs_block[i];
        for (0..count) |_|
            try blocks.append(.{
                .total_count = rs_block[i + 1],
                .data_count = rs_block[i + 2],
            });
    }
    return try blocks.toOwnedSlice();
}

test {
    const blocks = try rs_blocks(tst.allocator, 21, .L);
    defer tst.allocator.free(blocks);
    // std.debug.print("{any}\n", .{blocks});
}

pub const Polynomial = struct {
    num: []u8,
    allocator: Allocator,

    const Self = @This();

    pub fn init(allocator: Allocator, num: []const u8, shift: usize) !Self {
        // Handle empty array case
        if (num.len == 0) {
            return error.EmptyPolynomial;
        }

        // Find first non-zero coefficient
        var offset: usize = 0;
        while (offset < num.len and num[offset] == 0) {
            offset += 1;
        }

        // If all coefficients are zero, keep at least one zero
        if (offset == num.len) {
            offset = num.len - 1;
        }

        // Calculate final length
        const trimmed_len = num.len - offset;
        const final_len = trimmed_len + shift;

        // Allocate memory for coefficients
        const coeffs = try allocator.alloc(u8, final_len);

        // Copy non-zero coefficients
        @memcpy(coeffs[0..trimmed_len], num[offset..]);

        // Add shift zeros
        if (shift > 0) {
            @memset(coeffs[trimmed_len..], 0);
        }

        return Self{
            .num = coeffs,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.num);
    }

    // Indexing equivalent to __getitem__
    pub fn get(self: Self, index: usize) u8 {
        std.debug.assert(index < self.num.len);
        return self.num[index];
    }

    // Iterator support - returns slice for iteration
    pub fn iter(self: Self) []const u8 {
        return self.num;
    }

    // Length equivalent to __len__
    pub fn len(self: Self) usize {
        return self.num.len;
    }

    // Multiplication equivalent to __mul__
    pub fn mul(
        self: Self,
        other: Self,
    ) !Self {
        const result_len = self.len() + other.len() - 1;
        const result_coeffs = try self.allocator.alloc(u8, result_len);

        // Initialize result to zeros
        @memset(result_coeffs, 0);

        // Polynomial multiplication
        for (self.iter(), 0..) |self_coeff, i| {
            for (other.iter(), 0..) |other_coeff, j| {
                if (self_coeff != 0 and other_coeff != 0) {
                    result_coeffs[i + j] ^= gexp(glog(self_coeff) + glog(other_coeff));
                }
            }
        }

        return Self{
            .num = result_coeffs,
            .allocator = self.allocator,
        };
    }

    // Modulo equivalent to __mod__
    pub fn mod(self: Self, other: Self) !Self {
        const difference = @as(i32, @intCast(self.len())) - @as(i32, @intCast(other.len()));
        std.debug.print("{}\n", .{difference});
        if (difference < 0) {
            return error.Hi;
            // Return copy of self
            // const result_coeffs = try self.allocator.alloc(u8, self.len());
            // @memcpy(result_coeffs, self.num);
            // return Self{
            //     .num = result_coeffs,
            //     .allocator = self.allocator,
            // };
        }

        const ratio = @abs(@as(isize, glog(self.get(0))) - glog(other.get(0)));
        _ = ratio; // autofix
        const min_len = @min(self.len(), other.len());
        _ = min_len; // autofix
        const result_len = self.len();

        const result_coeffs = try self.allocator.alloc(u8, result_len);

        // // XOR operation for overlapping coefficients
        // for (0..min_len) |i| {
        //     result_coeffs[i] = self.get(i) ^ gexp(glog(other.get(i)) + ratio);
        // }

        // // Copy remaining coefficients if self is longer
        // if (difference > 0) {
        //     const remaining_start = other.len();
        //     @memcpy(result_coeffs[remaining_start..], self.num[remaining_start..]);
        // }

        var temp_poly = Self{
            .num = result_coeffs,
            .allocator = self.allocator,
        };
        // defer temp_poly.deinit();

        // Recursive call
        return temp_poly.mod(other);
    }

    // Helper method for debugging
    pub fn print(self: Self) void {
        std.debug.print("Polynomial coefficients: [", .{});
        for (self.num, 0..) |coeff, i| {
            if (i > 0) std.debug.print(", ", .{});
            std.debug.print("{}", .{coeff});
        }
        std.debug.print("]\n", .{});
    }
};

test "polynomial operations" {
    const allocator = std.testing.allocator;

    // Create polynomials
    var poly1 = try Polynomial.init(allocator, &[_]u8{ 1, 2, 3 }, 1);
    defer poly1.deinit();

    var poly2 = try Polynomial.init(allocator, &[_]u8{ 2, 1 }, 0);
    defer poly2.deinit();

    // Multiplication
    var result = try poly1.mul(poly2);
    defer result.deinit();

    // Modulo operation
    var mod_result = try poly1.mod(poly2);
    defer mod_result.deinit();

    // Iteration
    for (result.iter()) |coeff| {
        std.debug.print("{} ", .{coeff});
    }
}

const RSBlockTable = [_][]const u8{
    // L
    // M
    // Q
    // H
    // 1
    &.{ 1, 26, 19 },
    &.{ 1, 26, 16 },
    &.{ 1, 26, 13 },
    &.{ 1, 26, 9 },
    // 2
    &.{ 1, 44, 34 },
    &.{ 1, 44, 28 },
    &.{ 1, 44, 22 },
    &.{ 1, 44, 16 },
    // 3
    &.{ 1, 70, 55 },
    &.{ 1, 70, 44 },
    &.{ 2, 35, 17 },
    &.{ 2, 35, 13 },
    // 4
    &.{ 1, 100, 80 },
    &.{ 2, 50, 32 },
    &.{ 2, 50, 24 },
    &.{ 4, 25, 9 },
    // 5
    &.{ 1, 134, 108 },
    &.{ 2, 67, 43 },
    &.{ 2, 33, 15, 2, 34, 16 },
    &.{ 2, 33, 11, 2, 34, 12 },
    // 6
    &.{ 2, 86, 68 },
    &.{ 4, 43, 27 },
    &.{ 4, 43, 19 },
    &.{ 4, 43, 15 },
    // 7
    &.{ 2, 98, 78 },
    &.{ 4, 49, 31 },
    &.{ 2, 32, 14, 4, 33, 15 },
    &.{ 4, 39, 13, 1, 40, 14 },
    // 8
    &.{ 2, 121, 97 },
    &.{ 2, 60, 38, 2, 61, 39 },
    &.{ 4, 40, 18, 2, 41, 19 },
    &.{ 4, 40, 14, 2, 41, 15 },
    // 9
    &.{ 2, 146, 116 },
    &.{ 3, 58, 36, 2, 59, 37 },
    &.{ 4, 36, 16, 4, 37, 17 },
    &.{ 4, 36, 12, 4, 37, 13 },
    &.{ 6, 58, 36, 2, 59, 37 },
    &.{ 4, 46, 20, 6, 47, 21 },
    &.{ 7, 42, 14, 4, 43, 15 },
    // 13
    &.{ 4, 133, 107 },
    &.{ 8, 59, 37, 1, 60, 38 },
    &.{ 8, 44, 20, 4, 45, 21 },
    &.{ 12, 33, 11, 4, 34, 12 },
    // 14
    &.{ 3, 145, 115, 1, 146, 116 },
    &.{ 4, 64, 40, 5, 65, 41 },
    &.{ 11, 36, 16, 5, 37, 17 },
    &.{ 11, 36, 12, 5, 37, 13 },
    // 15
    &.{ 5, 109, 87, 1, 110, 88 },
    &.{ 5, 65, 41, 5, 66, 42 },
    &.{ 5, 54, 24, 7, 55, 25 },
    &.{ 11, 36, 12, 7, 37, 13 },
    // 16
    &.{ 5, 122, 98, 1, 123, 99 },
    &.{ 7, 73, 45, 3, 74, 46 },
    &.{ 15, 43, 19, 2, 44, 20 },
    &.{ 3, 45, 15, 13, 46, 16 },
    // 17
    &.{ 1, 135, 107, 5, 136, 108 },
    &.{ 10, 74, 46, 1, 75, 47 },
    &.{ 1, 50, 22, 15, 51, 23 },
    &.{ 2, 42, 14, 17, 43, 15 },
    // 18
    &.{ 5, 150, 120, 1, 151, 121 },
    &.{ 9, 69, 43, 4, 70, 44 },
    &.{ 17, 50, 22, 1, 51, 23 },
    &.{ 2, 42, 14, 19, 43, 15 },
    // 19
    &.{ 3, 141, 113, 4, 142, 114 },
    &.{ 3, 70, 44, 11, 71, 45 },
    &.{ 17, 47, 21, 4, 48, 22 },
    &.{ 9, 39, 13, 16, 40, 14 },
    // 20
    &.{ 3, 135, 107, 5, 136, 108 },
    &.{ 3, 67, 41, 13, 68, 42 },
    &.{ 15, 54, 24, 5, 55, 25 },
    &.{ 15, 43, 15, 10, 44, 16 },
    // 21
    &.{ 4, 144, 116, 4, 145, 117 },
    &.{ 17, 68, 42 },
    &.{ 17, 50, 22, 6, 51, 23 },
    &.{ 19, 46, 16, 6, 47, 17 },
    // 22
    &.{ 2, 139, 111, 7, 140, 112 },
    &.{ 17, 74, 46 },
    &.{ 7, 54, 24, 16, 55, 25 },
    &.{ 34, 37, 13 },
    // 23
    &.{ 4, 151, 121, 5, 152, 122 },
    &.{ 4, 75, 47, 14, 76, 48 },
    &.{ 11, 54, 24, 14, 55, 25 },
    &.{ 16, 45, 15, 14, 46, 16 },
    // 24
    &.{ 6, 147, 117, 4, 148, 118 },
    &.{ 6, 73, 45, 14, 74, 46 },
    &.{ 11, 54, 24, 16, 55, 25 },
    &.{ 30, 46, 16, 2, 47, 17 },
    // 25
    &.{ 8, 132, 106, 4, 133, 107 },
    &.{ 8, 75, 47, 13, 76, 48 },
    &.{ 7, 54, 24, 22, 55, 25 },
    &.{ 22, 45, 15, 13, 46, 16 },
    // 26
    &.{ 10, 142, 114, 2, 143, 115 },
    &.{ 19, 74, 46, 4, 75, 47 },
    &.{ 28, 50, 22, 6, 51, 23 },
    &.{ 33, 46, 16, 4, 47, 17 },
    // 27
    &.{ 8, 152, 122, 4, 153, 123 },
    &.{ 22, 73, 45, 3, 74, 46 },
    &.{ 8, 53, 23, 26, 54, 24 },
    &.{ 12, 45, 15, 28, 46, 16 },
    // 28
    &.{ 3, 147, 117, 10, 148, 118 },
    &.{ 3, 73, 45, 23, 74, 46 },
    &.{ 4, 54, 24, 31, 55, 25 },
    &.{ 11, 45, 15, 31, 46, 16 },
    // 29
    &.{ 7, 146, 116, 7, 147, 117 },
    &.{ 21, 73, 45, 7, 74, 46 },
    &.{ 1, 53, 23, 37, 54, 24 },
    &.{ 19, 45, 15, 26, 46, 16 },
    // 30
    &.{ 5, 145, 115, 10, 146, 116 },
    &.{ 19, 75, 47, 10, 76, 48 },
    &.{ 15, 54, 24, 25, 55, 25 },
    &.{ 23, 45, 15, 25, 46, 16 },
    // 31
    &.{ 13, 145, 115, 3, 146, 116 },
    &.{ 2, 74, 46, 29, 75, 47 },
    &.{ 42, 54, 24, 1, 55, 25 },
    &.{ 23, 45, 15, 28, 46, 16 },
    // 32
    &.{ 17, 145, 115 },
    &.{ 10, 74, 46, 23, 75, 47 },
    &.{ 10, 54, 24, 35, 55, 25 },
    &.{ 19, 45, 15, 35, 46, 16 },
    // 33
    &.{ 17, 145, 115, 1, 146, 116 },
    &.{ 14, 74, 46, 21, 75, 47 },
    &.{ 29, 54, 24, 19, 55, 25 },
    &.{ 11, 45, 15, 46, 46, 16 },
    // 34
    &.{ 13, 145, 115, 6, 146, 116 },
    &.{ 14, 74, 46, 23, 75, 47 },
    &.{ 44, 54, 24, 7, 55, 25 },
    &.{ 59, 46, 16, 1, 47, 17 },
    // 35
    &.{ 12, 151, 121, 7, 152, 122 },
    &.{ 12, 75, 47, 26, 76, 48 },
    &.{ 39, 54, 24, 14, 55, 25 },
    &.{ 22, 45, 15, 41, 46, 16 },
    // 36
    &.{ 6, 151, 121, 14, 152, 122 },
    &.{ 6, 75, 47, 34, 76, 48 },
    &.{ 46, 54, 24, 10, 55, 25 },
    &.{ 2, 45, 15, 64, 46, 16 },
    // 37
    &.{ 17, 152, 122, 4, 153, 123 },
    &.{ 29, 74, 46, 14, 75, 47 },
    &.{ 49, 54, 24, 10, 55, 25 },
    &.{ 24, 45, 15, 46, 46, 16 },
    // 38
    &.{ 4, 152, 122, 18, 153, 123 },
    &.{ 13, 74, 46, 32, 75, 47 },
    &.{ 48, 54, 24, 14, 55, 25 },
    &.{ 42, 45, 15, 32, 46, 16 },
    // 39
    &.{ 20, 147, 117, 4, 148, 118 },
    &.{ 40, 75, 47, 7, 76, 48 },
    &.{ 43, 54, 24, 22, 55, 25 },
    &.{ 10, 45, 15, 67, 46, 16 },
    // 40
    &.{ 19, 148, 118, 6, 149, 119 },
    &.{ 18, 75, 47, 31, 76, 48 },
    &.{ 34, 54, 24, 34, 55, 25 },
    &.{ 20, 45, 15, 61, 46, 16 },
};
