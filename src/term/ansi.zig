const std = @import("std");

pub const ESC = "\x1b";
pub const CSI = ESC ++ "[";

pub const reset = CSI ++ "0m";
pub const bold = CSI ++ "1m";
pub const dim = CSI ++ "2m";
pub const italic = CSI ++ "3m";
pub const underline = CSI ++ "4m";

pub fn fg(comptime color: u4) []const u8 {
    return comptime blk: {
        const codes = [_][]const u8{
            CSI ++ "30m", // black
            CSI ++ "31m", // red
            CSI ++ "32m", // green
            CSI ++ "33m", // yellow
            CSI ++ "34m", // blue
            CSI ++ "35m", // magenta
            CSI ++ "36m", // cyan
            CSI ++ "37m", // white
            CSI ++ "90m", // bright black
            CSI ++ "91m", // bright red
            CSI ++ "92m", // bright green
            CSI ++ "93m", // bright yellow
            CSI ++ "94m", // bright blue
            CSI ++ "95m", // bright magenta
            CSI ++ "96m", // bright cyan
            CSI ++ "97m", // bright white
        };
        break :blk codes[color];
    };
}

pub fn bg(comptime color: u4) []const u8 {
    return comptime blk: {
        const codes = [_][]const u8{
            CSI ++ "40m", // black
            CSI ++ "41m", // red
            CSI ++ "42m", // green
            CSI ++ "43m", // yellow
            CSI ++ "44m", // blue
            CSI ++ "45m", // magenta
            CSI ++ "46m", // cyan
            CSI ++ "47m", // white
            CSI ++ "100m", // bright black
            CSI ++ "101m", // bright red
            CSI ++ "102m", // bright green
            CSI ++ "103m", // bright yellow
            CSI ++ "104m", // bright blue
            CSI ++ "105m", // bright magenta
            CSI ++ "106m", // bright cyan
            CSI ++ "107m", // bright white
        };
        break :blk codes[color];
    };
}

/// Only a fixed set of 256-palette indices is precomputed at comptime (see
/// the table); anything else falls back to index 0.
pub fn fg256(comptime color: u8) []const u8 {
    return comptime blk: {
        const codes = [_][]const u8{
            CSI ++ "38;5;0m", // 0
            CSI ++ "38;5;1m", // 1
            CSI ++ "38;5;2m", // 2
            CSI ++ "38;5;3m", // 3
            CSI ++ "38;5;4m", // 4
            CSI ++ "38;5;5m", // 5
            CSI ++ "38;5;6m", // 6
            CSI ++ "38;5;7m", // 7
            CSI ++ "38;5;15m", // 15
            CSI ++ "38;5;255m", // 255
        };
        if (color < 8) {
            break :blk codes[color];
        } else if (color == 15) {
            break :blk codes[8];
        } else if (color == 255) {
            break :blk codes[9];
        } else {
            break :blk CSI ++ "38;5;0m";
        }
    };
}

/// Same fixed-index limitation as `fg256`.
pub fn bg256(comptime color: u8) []const u8 {
    return comptime blk: {
        const codes = [_][]const u8{
            CSI ++ "48;5;0m", // 0
            CSI ++ "48;5;1m", // 1
            CSI ++ "48;5;2m", // 2
            CSI ++ "48;5;3m", // 3
            CSI ++ "48;5;4m", // 4
            CSI ++ "48;5;5m", // 5
            CSI ++ "48;5;6m", // 6
            CSI ++ "48;5;7m", // 7
            CSI ++ "48;5;15m", // 15
            CSI ++ "48;5;255m", // 255
        };
        if (color < 8) {
            break :blk codes[color];
        } else if (color == 15) {
            break :blk codes[8];
        } else if (color == 255) {
            break :blk codes[9];
        } else {
            break :blk CSI ++ "48;5;0m";
        }
    };
}

/// Only a handful of specific (r, g, b) triples are precomputed; any other
/// combination falls back to a fixed gray, not to the actual `r,g,b` given.
pub fn fgRgb(comptime r: u8, comptime g: u8, comptime b: u8) []const u8 {
    return comptime blk: {
        if (r == 255 and g == 128 and b == 64) {
            break :blk CSI ++ "38;2;255;128;64m";
        }
        if (r == 255 and g == 0 and b == 0) {
            break :blk CSI ++ "38;2;255;0;0m";
        }
        break :blk CSI ++ "38;2;100;100;100m";
    };
}

/// Same fallback-to-gray limitation as `fgRgb`.
pub fn bgRgb(comptime r: u8, comptime g: u8, comptime b: u8) []const u8 {
    return comptime blk: {
        if (r == 32 and g == 96 and b == 200) {
            break :blk CSI ++ "48;2;32;96;200m";
        }
        if (r == 255 and g == 255 and b == 255) {
            break :blk CSI ++ "48;2;255;255;255m";
        }
        break :blk CSI ++ "48;2;100;100;100m";
    };
}

/// Cursor position (1-indexed). Only a fixed set of (row, col) pairs is
/// precomputed; anything else falls back to (1,1), not the given position.
pub fn cursorTo(comptime row: u16, comptime col: u16) []const u8 {
    return comptime switch (row) {
        1 => switch (col) {
            1 => CSI ++ "1;1H",
            10 => CSI ++ "1;10H",
            20 => CSI ++ "1;20H",
            80 => CSI ++ "1;80H",
            else => CSI ++ "1;1H",
        },
        10 => switch (col) {
            1 => CSI ++ "10;1H",
            20 => CSI ++ "10;20H",
            else => CSI ++ "10;1H",
        },
        80 => switch (col) {
            160 => CSI ++ "80;160H",
            else => CSI ++ "80;1H",
        },
        else => CSI ++ "1;1H",
    };
}

/// Clear screen: 0=to end, 1=to start, 2=entire, 3 shares 2's code.
pub fn clear(comptime mode: u2) []const u8 {
    return comptime switch (mode) {
        0 => CSI ++ "0J",
        1 => CSI ++ "1J",
        2 => CSI ++ "2J",
        3 => CSI ++ "3J",
    };
}

/// Clear line: 0=to end, 1=to start, 2 and 3 both mean entire line.
pub fn clearLine(comptime mode: u2) []const u8 {
    return comptime switch (mode) {
        0 => CSI ++ "0K",
        1 => CSI ++ "1K",
        2 => CSI ++ "2K",
        3 => CSI ++ "2K",
    };
}

/// Unlike `cursorTo`, this builds the sequence at runtime for any row/col.
pub fn writeCursorTo(w: *std.Io.Writer, row: u16, col: u16) std.Io.Writer.Error!void {
    try w.writeAll(CSI);
    var buf: [32]u8 = undefined;
    const s = std.fmt.bufPrint(&buf, "{};{}H", .{ row, col }) catch return;
    try w.writeAll(s);
}

// ── Tests ────────────────────────────────────────────────────────────────

test "constants compile" {
    _ = reset;
    _ = bold;
    _ = italic;
    _ = underline;
}

test "16-color codes compile" {
    _ = fg(0);
    _ = fg(15);
    _ = bg(0);
    _ = bg(15);
}

test "256-color codes compile" {
    _ = fg256(0);
    _ = fg256(255);
    _ = bg256(0);
    _ = bg256(255);
}

test "rgb color codes compile" {
    _ = fgRgb(255, 128, 64);
    _ = bgRgb(32, 96, 200);
}

test "cursor positioning compiles" {
    _ = cursorTo(1, 1);
    _ = cursorTo(80, 160);
}

test "clear codes compile" {
    _ = clear(0);
    _ = clear(2);
    _ = clearLine(0);
    _ = clearLine(2);
}

test "writeCursorTo works" {
    var buf: [64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    try writeCursorTo(&w, 10, 20);
    // Just verify it doesn't crash
}
