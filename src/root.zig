const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;
pub const Astro = @import("astro.zig");

const koino = @import("koino");

test {
    std.testing.refAllDecls(@This());
}
test {
    var k = try koino.parser.Parser.init(tst.allocator, .{
        .extensions = .{},
    });
    defer k.deinit();

    try k.feed("**Hello**, [world](http://world.gov)!");

    const doc = try k.finish();
    var iter = doc.traverseIterator();
    while (iter.next()) |node| {
        switch (node) {
            else => {
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

    // std.debug.print("{s}\n", .{output});
}
