const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;

const Zyrtex = @import("zyrtex/ast.zig");

pub const koino = @import("koino");

test {
    std.testing.refAllDecls(@This());

    var k = try koino.parser.Parser.init(tst.allocator, koino.Options{
        .extensions = .{
            .autolink = true,
            .strikethrough = true,
            .table = true,
            .tagfilter = true,
        },
    });
    defer k.deinit();

    try k.feed(
        \\## heading 1
        \\## heading 2
        \\```mermaid
        \\graph TD
        \\A-->B
        \\```
    );

    const doc = try k.finish();
    var iter = doc.traverseIterator();
    while (iter.next()) |node| {
        switch (node) {
            .Start => {
                std.debug.print("Node: {}\n", .{node.Start.data.value});
            },
            .End => {
                // std.debug.print("Node: {}\n", .{n});
            },
        }
    }

    const output = blk: {
        var arr = std.ArrayList(u8).init(tst.allocator);
        errdefer arr.deinit();
        try koino.html.print(arr.writer(), tst.allocator, .{}, doc);
        break :blk try arr.toOwnedSlice();
    };
    defer tst.allocator.free(output);
    defer doc.deinit();

    std.debug.print("{s}\n", .{output});
}
