const std = @import("std");
const Image = @import("../image.zig");

pub const Header = extern struct {
    signature: u16,
    file_size: u32 = 0,
    reserved: u32 = 0,
    offset: u16,
};
comptime {
    // @compileLog(@sizeOf(Header));
    // std.debug.assert(@sizeOf(Header) == 14);
}

pub const InfoHeader = extern struct {
    info_size: u32 = 40,
    width: u32,
    height: u32,
    planes: u16 = 1,
    bits_per_pixel: u16 = 24,
    compression: u32 = 0,
    image_size: u32 = 0,
    x_pixels_per_m: u32,
    y_pixels_per_m: u32,
    colors_used: u32,
    important_colors: u32 = 0,
};

pub const ColorTable = packed struct {
    colors: [*]Color,
};

const Color = extern struct {
    r: u8,
    g: u8,
    b: u8,
    _: u8,
};

pub const File = extern struct {
    header: Header,
    info: InfoHeader,
    colors: ColorTable,
    pixels: [*]Color,
};

pub fn read(_: std.mem.Allocator, path: []const u8) !Image {
    var file = try std.fs.cwd().openFile(path, .{
        .mode = .read_only,
    });
    defer file.close();

    // const contents = try file.readToEndAlloc(alloc, 100_000_000);
    // var reader = std.io.bufferedReader(file.reader());
    var reader = file.reader();

    const header = try reader.readStruct(Header);
    std.debug.print("{} {any}\n", .{ header.signature, header });
    const info = try reader.readStruct(InfoHeader);
    std.debug.print("{any}\n", .{info});

    return Image{};
}

pub fn write(_: Image, _: []const u8) !void {}

test {
    // _ = try read(std.testing.allocator, "img.bmp");
}
