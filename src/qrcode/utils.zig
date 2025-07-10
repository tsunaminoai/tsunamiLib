const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;
const code_point = @import("code_point");
const ascii = std.ascii;
pub const Mode = enum {
    number,
    alpha,
    byte,
    kanji,
};
pub const Penalty = enum(u8) {
    n1 = 3,
    n2 = 3,
    n3 = 40,
    n4 = 10,
    pub const Info = struct {
        horizRuns: Array(LinearRun),
        vertRuns: Array(LinearRun),
        twoByTwoBoxes: Array(@TypeOf(.{ usize, usize })),
        horizFalseFinders: Array(LinearRun),
        vertFalseFinders: Array(LinearRun),
        numDarkModules: usize = 0,
        penaltyPoints: Array(@TypeOf(.{ usize, usize, usize, usize })),
    };
};
pub const LinearRun = struct {
    startX: usize = 0,
    startY: usize = 0,
    runLen: usize = 0,
};

pub const Segment = extern struct {
    mode: u4, // The segment mode is always a 4-bit field.
    count: u8, // The character count’s field width depends on the mode and version.
    data: [17]u8,
    terminator: u4 = 0, // The terminator is normally four “0” bits, but fewer if the data codeword capacity is reached.
    //The bit padding is between zero to seven “0” bits, to fill all unused bits in the last byte.
    // The byte padding consists of alternating (hexadecimal) EC and 11 until the capacity is reached.
};
pub const Block = struct {};
pub const ECCLevel = union(enum) {
    low,
    med,
    qrt,
    hi,
};
pub const Codeword = struct {
    value: u8,
    // preInterleaveIndex
    // blockIndex
    // indexInBlock
    // postInterleaveIndex
    pub fn init(v: u8) !Codeword {
        if (v > 255) return error.InvalidValue;
        return .{
            .value = v,
        };
    }
};

pub const ECCCodeWordsPerBlock = [_][]const i16{
    &.{ -1, 7, 10, 15, 20, 26, 18, 20, 24, 30, 18, 20, 24, 26, 30, 22, 24, 28, 30, 28, 28, 28, 28, 30, 30, 26, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 },
    &.{ -1, 10, 16, 26, 18, 24, 16, 18, 22, 22, 26, 30, 22, 22, 24, 24, 28, 28, 26, 26, 26, 26, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28 },
    &.{ -1, 13, 22, 18, 26, 18, 24, 18, 22, 20, 24, 28, 26, 24, 20, 30, 24, 28, 28, 26, 30, 28, 30, 30, 30, 30, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 },
    &.{ -1, 17, 28, 22, 16, 22, 28, 26, 26, 24, 28, 24, 28, 22, 24, 24, 30, 28, 28, 26, 28, 30, 24, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 },
};
pub const NumECCBlocks = [_][]const i16{
    &.{ -1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 7, 8, 8, 9, 9, 10, 12, 12, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21, 22, 24, 25 },
    &.{ -1, 1, 1, 1, 2, 2, 4, 4, 4, 5, 5, 5, 8, 9, 9, 10, 10, 11, 13, 14, 16, 17, 17, 18, 20, 21, 23, 25, 26, 28, 29, 31, 33, 35, 37, 38, 40, 43, 45, 47, 49 },
    &.{ -1, 1, 1, 2, 2, 4, 4, 6, 6, 8, 8, 8, 10, 12, 16, 12, 17, 16, 18, 21, 20, 23, 23, 25, 27, 29, 34, 34, 35, 38, 40, 43, 45, 48, 51, 53, 56, 59, 62, 65, 68 },
    &.{ -1, 1, 1, 2, 4, 4, 4, 5, 6, 8, 8, 11, 11, 16, 16, 18, 16, 19, 21, 25, 25, 25, 34, 30, 32, 35, 37, 40, 42, 45, 48, 51, 54, 57, 60, 63, 66, 70, 74, 77, 81 },
};

pub const Version = enum(u8) {
    invalid = 0,
    @"1" = 1,
    @"40" = 40,
    pub fn get_config(self: Version) Config {
        return switch (@intFromEnum(self)) {
            1...9 => .{ .bits = 148, .codewords = 20 },
            10...26 => .{ .bits = 156, .codewords = 20 },
            27...40 => .{ .bits = 156, .codewords = 20 },
            else => unreachable,
        };
    }
    pub const Config = struct {
        ecc: enum { L, M, Q, H, check } = .check,
        bits: usize,
        codewords: usize,
    };
};

pub fn analyze_unicode(input: []const u8) !Mode {
    var tests = packed struct(u4) {
        is_num: bool = false,
        is_alpha: bool = false,
        is_byte: bool = false,
        is_kanji: bool = false,
    }{};
    var iter: code_point.Iterator = .init(input);
    while (iter.next()) |cp| {
        // if (!uc.utf8ValidCodepoint(cp)) return error.InvalidUnicode;
        if (cp.code < 255) {
            const ch: u8 = @intCast(cp.code);
            if (ascii.isAlphanumeric(ch)) {
                if (ascii.isAlphabetic(ch))
                    tests.is_alpha = true
                else
                    tests.is_num = true;
            } else tests.is_byte = true;
        }
    }
    return switch (@as(u4, @bitCast(tests))) {
        0b0001 => .number,
        0b0010...0b0011 => .alpha,
        0b0100...0b0111 => .byte,
        0b1000...0b1111 => .kanji,
        else => error.InvalidUnicode,
    };
}

pub const Module = union(enum) {
    unfilled: Unfilled,
    function: Function,
    code_word: CodeWord,
    remainder: Remainder,
    mask: Mask,
    filled: Filled,

    pub const Filled = struct {
        is_new: bool,
    };
    pub const Unfilled = struct {};
    pub const Function = struct {
        pub const Kind = enum {
            finder,
            separator,
            timing,
            alignment,
            format_info,
            version_info,
            dark,
        };
    };
    pub const CodeWord = struct {};
    pub const Remainder = struct {};
    pub const Mask = struct {};
};

pub const QRCode = struct {
    side_len: usize,
    modules: Array(Array(Module)), // side x side
    version: Version,
    ecc_level: ECCLevel,

    allocator: Allocator,

    pub fn clearNewFlags(self: *QRCode) void {
        for (self.modules.items) |column| {
            for (column.items) |mod| {
                switch (mod) {
                    .filled => |*f| f.is_new = false,
                    else => continue,
                }
            }
        }
    }

    pub fn makeZigZagScan(self: QRCode) ![]@TypeOf(.{ usize, usize }) {
        var res = Array(@TypeOf(.{ usize, usize })).init(self.allocator);
        defer res.deinit();
        var right: isize = self.side_len;
        while (right >= 1) : (right -= 2) {
            if (right == 6) right = 5;
            var vert: usize = 0;
            while (vert < self.side_len) : (vert += 1) {
                for (0..2) |j| {
                    const x = right - j;
                    const up = ((right + 1) & 2) == 0;
                    const y = if (up) self.side_len - 1 - vert else vert;
                    if (self.modules.items[x].items[y] == .function)
                        try res.append(.{ x, y });
                }
            }
        }
        return try res.toOwnedSlice();
    }
};

// from https://www.nayuki.io/page/creating-a-qr-code-step-by-step
test "hello world" {
    const input = "Hello, world! 123";
    try tst.expectEqual(.byte, try analyze_unicode(input));
    const cfg = Version.@"1".get_config();
    std.debug.print("{}\n", .{cfg});
    const expect_codewords = [_]u8{ 0x41, 0x14, 0x86, 0x56, 0xC6, 0xC6, 0xF2, 0xC2, 0x07, 0x76, 0xF7, 0x26, 0xC6, 0x42, 0x12, 0x03, 0x13, 0x23, 0x30, 0x85, 0xA9, 0x5E, 0x07, 0x0A, 0x36, 0xC9 };
    _ = expect_codewords; // autofix
}
