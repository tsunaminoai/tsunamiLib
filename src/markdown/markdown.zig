const std = @import("std");
const Allocator = std.mem.Allocator;
const Array = std.ArrayList;
/// Heres what I'd like to see:
/// MD -> Lines -> AST
/// Needs to have formatting for each line type
/// Assuming no multiline elements for now
fn eatWhitespace(txt: []const u8, idx: usize) usize {
    for (txt[idx..], 0..) |c, i| {
        if (!std.ascii.isWhitespace(c)) {
            return i + idx;
        }
    }
}

const LineType = enum(u8) {
    heading = '#',
    block_quote = '>',
    list_u = '-',
    text,
    blank,

    fn toType(c: u8) LineType {
        return std.meta.intToEnum(LineType, c) catch return .text;
    }
};

pub const MDDoc = struct {
    lines: []Line,
    const Line = struct {
        txt: []const u8,
        type: LineType,
        pub fn parse(src: []const u8) !Line {
            return .{
                .txt = src,
                .type = if (src.len > 0) LineType.toType(src[0]) else .blank,
            };
        }
    };
    fn parse(allocator: std.mem.Allocator, md: []const u8) !MDDoc {
        var idx: usize = 0;
        var lines = Array(MDDoc.Line).init(allocator);
        var split = std.mem.splitAny(u8, md, "\n");
        while (split.next()) |l| {
            const line = try Line.parse(l);
            try lines.append(line);
            idx += line.txt.len;
        }
        return MDDoc{
            .lines = try lines.toOwnedSlice(),
        };
    }
    pub fn format(
        self: MDDoc,
        _: []const u8,
        _: anytype,
        writer: anytype,
    ) !void {
        for (self.lines) |l| {
            try writer.print("{} {s}\n", .{ l.type, l.txt });
        }
    }
};

const testDoc =
    \\# H1
    \\
    \\Paragaph
    \\More text in the same paragraph
    \\
    \\New Paragraph
    \\
    \\> Quote
;

const expecteds = .{
    LineType.heading,
    LineType.blank,
    LineType.text,
    LineType.text,
    LineType.blank,
    LineType.text,
    LineType.blank,
    LineType.block_quote,
};
const testing = std.testing;
const expectEql = std.testing.expectEqualDeep;
test "MDDoc" {
    const md = try MDDoc.parse(testing.allocator, testDoc[0..]);
    defer testing.allocator.free(md.lines);

    inline for (expecteds, 0..) |exp, i| {
        try expectEql(md.lines[i].type, exp);
    }
}
