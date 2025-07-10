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

// from https://www.nayuki.io/page/creating-a-qr-code-step-by-step
test "hello world" {
    const input = "Hello, world! 123";
    try tst.expectEqual(.byte, try analyze_unicode(input));
}
