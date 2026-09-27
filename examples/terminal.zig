//! `term.ansi` comptime colored strings, `term.braille` rendering a small
//! sparkline, and `term.size.get` reporting the terminal size (or "not a
//! tty", the expected case when stdout is a pipe under `zig build`).
const std = @import("std");
const ts = @import("tsunami");
const term = ts.term;

/// One braille cell packs a 2x4 dot grid; here each column of the sparkline
/// uses only the left dot column, so 4 sample levels map to 4 dot rows.
fn levelDots(level: u2) term.braille.BrailleSet {
    return switch (level) {
        0 => .blank,
        1 => term.braille.BrailleSet.bottom_left,
        2 => term.braille.BrailleSet.bottom_left.mix(.middle_filled),
        3 => .full,
    };
}

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    try w.print("{s}green{s} {s}bold{s}\n", .{ term.ansi.fg(2), term.ansi.reset, term.ansi.bold, term.ansi.reset });

    const samples = [_]u8{ 1, 3, 6, 9, 5, 2, 8, 9 };
    var max: u8 = 0;
    for (samples) |s| max = @max(max, s);
    var line: [samples.len * 3]u8 = undefined;
    var pos: usize = 0;
    for (samples) |s| {
        const level: u2 = @intCast(@min(3, s * 4 / (max + 1)));
        const cp = levelDots(level).asUnicode();
        const n = std.unicode.utf8Encode(cp, line[pos..]) catch unreachable;
        pos += n;
    }
    try w.print("sparkline: {s}\n", .{line[0..pos]});
    if (pos == 0) return error.ExampleFailed;

    if (term.size.get(std.posix.STDOUT_FILENO)) |size| {
        try w.print("terminal size: {d}x{d}\n", .{ size.cols, size.rows });
    } else {
        try w.print("terminal size: not a tty\n", .{});
    }
    try w.flush();
}
