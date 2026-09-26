const std = @import("std");

pub const contract = @import("contract.zig");

pub const dsp = struct {
    pub const fft = @import("dsp/fft.zig");
};

test {
    refAllNamespaces(@This());
}

fn refAllNamespaces(comptime T: type) void {
    inline for (comptime std.meta.declarations(T)) |d| {
        const v = @field(T, d.name);
        if (@TypeOf(v) == type and @typeInfo(v) == .@"struct" and std.meta.declarations(v).len > 0 and v != std) {
            std.testing.refAllDecls(v);
        }
    }
}
