const std = @import("std");
pub const Bitmap = @import("image/bmp.zig");

pub const Color = union {
    pub const RGB = struct {
        r: u8 = 0,
        g: u8 = 0,
        b: u8 = 0,
    };
    pub const RGBA = struct {
        r: u8 = 0,
        g: u8 = 0,
        b: u8 = 0,
        a: u8 = 0,
    };

    rgb: RGB,
    rgba: RGBA,
};

test {
    std.testing.refAllDecls(@This());
}
