/// Braille set utilities for rendering tracks and visualizations.
/// Allocation-free — the caller supplies the byte buffer the UTF-8 output
/// is rendered into.
const std = @import("std");

pub const BrailleSet = enum(u8) {
    blank = 0x00,
    one = 0x01,
    two = 0x02,
    three = 0x04,
    four = 0x08,
    five = 0x10,
    six = 0x20,
    seven = 0x40,
    eight = 0x80,
    full = 0xFF,
    _,

    pub const top_left: BrailleSet = .one;
    pub const top_right: BrailleSet = .four;
    pub const bottom_left: BrailleSet = .seven;
    pub const bottom_right: BrailleSet = .eight;
    pub const middle_filled: BrailleSet = @enumFromInt(@intFromEnum(BrailleSet.two) |
        @intFromEnum(BrailleSet.three) | @intFromEnum(BrailleSet.five) | @intFromEnum(BrailleSet.six));

    const base: u21 = 0x2800;

    pub fn mix(self: BrailleSet, other: BrailleSet) BrailleSet {
        return @enumFromInt(@intFromEnum(self) | @intFromEnum(other));
    }

    pub fn asUnicode(self: BrailleSet) u21 {
        return base + @as(u21, @intFromEnum(self));
    }
};

// ── Tests ────────────────────────────────────────────────────────────────

const tst = std.testing;

test "braille mix combines dots" {
    const combined = BrailleSet.top_left.mix(BrailleSet.bottom_right);
    try tst.expectEqual(@intFromEnum(BrailleSet.top_left) | @intFromEnum(BrailleSet.bottom_right), @intFromEnum(combined));
}

test "braille unicode encoding" {
    const blank_code = BrailleSet.blank.asUnicode();
    const full_code = BrailleSet.full.asUnicode();
    try tst.expectEqual(@as(u21, 0x2800), blank_code);
    try tst.expectEqual(@as(u21, 0x28FF), full_code);
}

test "braille middle_filled dot pattern" {
    const middle = BrailleSet.middle_filled;
    // Should be 0x02 | 0x04 | 0x10 | 0x20 = 0x36
    const expected = 0x02 | 0x04 | 0x10 | 0x20;
    try tst.expectEqual(expected, @intFromEnum(middle));
}
