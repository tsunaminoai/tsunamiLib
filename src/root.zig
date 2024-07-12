pub const termsize = @import("termsize.zig");

test {
    const std = @import("std");
    std.testing.refAllDecls(@This());
}
