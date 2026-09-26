//! Generic zero-copy helpers for reading/writing fixed-layout binary records
//! out of raw byte slices. `T` must be an `extern struct` so its Zig layout
//! matches its on-disk layout; endian conversion swaps every integer field.

const std = @import("std");
const contract = @import("../contract.zig");

pub const Error = error{TooSmall};

fn assertExternStruct(comptime T: type) void {
    const info = @typeInfo(T);
    if (info != .@"struct" or info.@"struct".layout != .@"extern")
        @compileError(@typeName(T) ++ ": blocks.view/read require an extern struct");
}

/// Zero-copy reinterpret of `bytes` as `*align(1) const T`. `T` must be an
/// `extern struct`. Errors (rather than panics) if `bytes` is too short.
pub fn view(comptime T: type, bytes: []const u8) Error!*align(1) const T {
    comptime assertExternStruct(T);
    if (bytes.len < @sizeOf(T)) return Error.TooSmall;
    return std.mem.bytesAsValue(T, bytes[0..@sizeOf(T)]);
}

/// Reads a `T` by value out of `bytes`, byte-swapping every integer field
/// (recursively, including nested structs and arrays) from `endian` to
/// native. Non-integer fields are handled the same way `std.mem.
/// byteSwapAllFields` treats them: enums swap their backing integer, floats
/// swap their raw bits, bools and zero-sized types pass through untouched.
pub fn read(comptime T: type, bytes: []const u8, comptime endian: std.builtin.Endian) Error!T {
    comptime assertExternStruct(T);
    if (bytes.len < @sizeOf(T)) return Error.TooSmall;
    var value = std.mem.bytesToValue(T, bytes[0..@sizeOf(T)]);
    if (comptime endian != @import("builtin").cpu.arch.endian())
        std.mem.byteSwapAllFields(T, &value);
    return value;
}

/// Inverse of `read`: byte-swaps `value` from native to `endian` (if needed)
/// and writes its raw bytes into `out`. Errors if `out` is too small.
pub fn writeInto(comptime T: type, value: T, comptime endian: std.builtin.Endian, out: []u8) Error!void {
    comptime assertExternStruct(T);
    if (out.len < @sizeOf(T)) return Error.TooSmall;
    var v = value;
    if (comptime endian != @import("builtin").cpu.arch.endian())
        std.mem.byteSwapAllFields(T, &v);
    @memcpy(out[0..@sizeOf(T)], std.mem.asBytes(&v));
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

const Inner = extern struct {
    a: u16,
    b: [3]u8,
};

const Rec = extern struct {
    id: u32,
    tag: u16,
    inner: Inner,
    counts: [2]u32,
};

fn sampleRec() Rec {
    return .{
        .id = 0xDEADBEEF,
        .tag = 0xABCD,
        .inner = .{ .a = 0x1122, .b = .{ 1, 2, 3 } },
        .counts = .{ 0x11223344, 0x55667788 },
    };
}

test "view rejects short slices and returns a zero-copy pointer" {
    const rec = sampleRec();
    const bytes = std.mem.asBytes(&rec);

    const v = try view(Rec, bytes);
    try testing.expectEqual(rec.id, v.id);
    try testing.expectEqual(rec.inner.a, v.inner.a);
    try testing.expectEqualSlices(u8, &rec.inner.b, &v.inner.b);

    try testing.expectError(Error.TooSmall, view(Rec, bytes[0 .. bytes.len - 1]));
}

test "read/writeInto round trip, native endian" {
    const rec = sampleRec();
    var buf: [@sizeOf(Rec)]u8 = undefined;
    try writeInto(Rec, rec, .little, &buf);
    const got = try read(Rec, &buf, .little);
    try testing.expectEqual(rec, got);
}

test "read/writeInto round trip, foreign endian swaps every integer field" {
    const rec = sampleRec();
    const foreign: std.builtin.Endian = comptime if (@import("builtin").cpu.arch.endian() == .little) .big else .little;

    var buf: [@sizeOf(Rec)]u8 = undefined;
    try writeInto(Rec, rec, foreign, &buf);
    try testing.expect(!std.mem.eql(u8, std.mem.asBytes(&rec), &buf));

    const got = try read(Rec, &buf, foreign);
    try testing.expectEqual(rec, got);
}

test "read errors on truncated input instead of panicking" {
    const rec = sampleRec();
    const bytes = std.mem.asBytes(&rec);
    try testing.expectError(Error.TooSmall, read(Rec, bytes[0 .. bytes.len - 1], .little));
}

test "writeInto errors on undersized output" {
    const rec = sampleRec();
    var buf: [@sizeOf(Rec) - 1]u8 = undefined;
    try testing.expectError(Error.TooSmall, writeInto(Rec, rec, .little, &buf));
}

test "fuzz read never panics on arbitrary bytes" {
    try std.testing.fuzz({}, struct {
        fn f(_: void, s: *std.testing.Smith) !void {
            var buf: [512]u8 = undefined;
            const n = s.slice(&buf);
            _ = read(Rec, buf[0..n], .little) catch return;
        }
    }.f, .{});
}
