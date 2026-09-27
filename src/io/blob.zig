//! A versioned, fixed-size save blob: `magic | version | payload | crc32`.
//! `Blob(T, magic, version)` derives T's byte layout at comptime by walking
//! `@typeInfo(T)`; the layout is little-endian for every integer/enum
//! regardless of target, so a blob written on one target decodes on another.
//! `decode` never trusts its input: a short buffer, wrong magic/version, or a
//! flipped bit anywhere in the payload comes back as a typed error rather
//! than a partially-applied `T` or an out-of-bounds read.
//!
//! Supported field types (recursive via `@typeInfo`): integers, bools, enums
//! (as their tag integer), fixed-size arrays, and nested structs. Anything
//! else (pointers, slices, unions, optionals, floats) is a `@compileError` at
//! `Blob` instantiation — a save format must have one unambiguous byte shape.

const std = @import("std");
const contract = @import("../contract.zig");

const header_len = 4 + 2; // magic + version
const trailer_len = 4; // crc32

/// A versioned fixed-size encoder/decoder for `T`.
///
/// `magic` identifies the blob kind; `version` is bumped for any layout
/// change that is not a pure encode/decode-compatible append. A mismatched
/// magic or version is rejected rather than guessed at.
pub fn Blob(comptime T: type, comptime magic: [4]u8, comptime version: u16) type {
    const payload_len = comptime encodedSize(T);
    return struct {
        /// Total encoded size: magic(4) + version(2) + payload + crc32(4).
        pub const size = header_len + payload_len + trailer_len;

        pub const Error = error{
            TooShort,
            BadMagic,
            BadVersion,
            BadChecksum,
        };

        /// Write `value` into `out` as `magic | version | payload | crc32`.
        pub fn encode(value: T, out: *[size]u8) void {
            out[0..4].* = magic;
            std.mem.writeInt(u16, out[4..6], version, .little);
            encodeValue(T, value, out[header_len .. header_len + payload_len]);
            const crc = std.hash.Crc32.hash(out[0 .. header_len + payload_len]);
            std.mem.writeInt(u32, out[header_len + payload_len ..][0..4], crc, .little);
        }

        /// Parse a blob, validating length, magic, version and checksum
        /// before decoding the payload. Never reads out of bounds and never
        /// partially applies a bad buffer.
        pub fn decode(bytes: []const u8) Error!T {
            if (bytes.len != size) return Error.TooShort;
            if (!std.mem.eql(u8, bytes[0..4], &magic)) return Error.BadMagic;
            const v = std.mem.readInt(u16, bytes[4..6], .little);
            if (v != version) return Error.BadVersion;
            const crc_expected = std.mem.readInt(u32, bytes[header_len + payload_len ..][0..4], .little);
            const crc_actual = std.hash.Crc32.hash(bytes[0 .. header_len + payload_len]);
            if (crc_actual != crc_expected) return Error.BadChecksum;
            return decodeValue(T, bytes[header_len .. header_len + payload_len]);
        }
    };
}

/// Comptime byte size of `T`'s encoded payload.
fn encodedSize(comptime T: type) usize {
    return switch (@typeInfo(T)) {
        .bool => 1,
        .int => |i| (i.bits + 7) / 8,
        .@"enum" => |e| encodedSize(e.tag_type),
        .array => |a| encodedSize(a.child) * a.len,
        .@"struct" => |s| blk: {
            var total: usize = 0;
            inline for (s.fields) |f| total += encodedSize(f.type);
            break :blk total;
        },
        else => @compileError("blob: unsupported field type " ++ @typeName(T) ++
            " (only bool, int, enum, fixed-size array, and nested struct are supported)"),
    };
}

fn encodeValue(comptime T: type, value: T, out: []u8) void {
    switch (@typeInfo(T)) {
        .bool => out[0] = @intFromBool(value),
        .int => |i| {
            const bytes = comptime encodedSize(T);
            const Same = std.meta.Int(.unsigned, i.bits);
            const U = std.meta.Int(.unsigned, bytes * 8);
            const widened: U = @as(Same, @bitCast(value));
            std.mem.writeInt(U, out[0..bytes], widened, .little);
        },
        .@"enum" => |e| encodeValue(e.tag_type, @intFromEnum(value), out),
        .array => |a| {
            const elem_size = encodedSize(a.child);
            for (value, 0..) |elem, idx| {
                encodeValue(a.child, elem, out[idx * elem_size .. (idx + 1) * elem_size]);
            }
        },
        .@"struct" => |s| {
            var off: usize = 0;
            inline for (s.fields) |f| {
                const fsize = encodedSize(f.type);
                encodeValue(f.type, @field(value, f.name), out[off .. off + fsize]);
                off += fsize;
            }
        },
        else => @compileError("blob: unsupported field type " ++ @typeName(T)),
    }
}

fn decodeValue(comptime T: type, in: []const u8) T {
    switch (@typeInfo(T)) {
        .bool => return in[0] != 0,
        .int => |i| {
            const bytes = comptime encodedSize(T);
            const Same = std.meta.Int(.unsigned, i.bits);
            const U = std.meta.Int(.unsigned, bytes * 8);
            const widened = std.mem.readInt(U, in[0..bytes], .little);
            const truncated: Same = @intCast(widened);
            return @bitCast(truncated);
        },
        .@"enum" => |e| {
            const tag = decodeValue(e.tag_type, in);
            return @enumFromInt(tag);
        },
        .array => |a| {
            var result: T = undefined;
            const elem_size = encodedSize(a.child);
            for (&result, 0..) |*elem, idx| {
                elem.* = decodeValue(a.child, in[idx * elem_size .. (idx + 1) * elem_size]);
            }
            return result;
        },
        .@"struct" => |s| {
            var result: T = undefined;
            var off: usize = 0;
            inline for (s.fields) |f| {
                const fsize = encodedSize(f.type);
                @field(result, f.name) = decodeValue(f.type, in[off .. off + fsize]);
                off += fsize;
            }
            return result;
        },
        else => @compileError("blob: unsupported field type " ++ @typeName(T)),
    }
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

const Color = enum(u8) { red, green, blue };

const Inner = struct {
    x: i16,
    y: i16,
};

const Sample = struct {
    id: u32,
    active: bool,
    color: Color,
    pos: Inner,
    tags: [3]u8,
};

const SampleBlob = Blob(Sample, .{ 'S', 'M', 'P', 'L' }, 1);

test "size accounts for header, payload and crc trailer" {
    // id(4) + active(1) + color(1) + pos(2+2) + tags(3) = 13 payload bytes.
    try testing.expectEqual(@as(usize, 4 + 2 + 13 + 4), SampleBlob.size);
}

test Blob {
    const value: Sample = .{
        .id = 0xDEADBEEF,
        .active = true,
        .color = .blue,
        .pos = .{ .x = -100, .y = 200 },
        .tags = .{ 1, 2, 3 },
    };
    var buf: [SampleBlob.size]u8 = undefined;
    SampleBlob.encode(value, &buf);
    const back = try SampleBlob.decode(&buf);
    try testing.expectEqual(value, back);
}

test "encode is deterministic and little-endian" {
    const value: Sample = .{
        .id = 1,
        .active = false,
        .color = .green,
        .pos = .{ .x = 0, .y = 0 },
        .tags = .{ 0, 0, 0 },
    };
    var buf: [SampleBlob.size]u8 = undefined;
    SampleBlob.encode(value, &buf);
    try testing.expectEqualSlices(u8, "SMPL", buf[0..4]);
    try testing.expectEqual(@as(u16, 1), std.mem.readInt(u16, buf[4..6], .little));
    // id = 1 little-endian at offset 6.
    try testing.expectEqual(@as(u32, 1), std.mem.readInt(u32, buf[6..10], .little));
}

test "decode rejects a too-short buffer" {
    var buf: [SampleBlob.size]u8 = undefined;
    SampleBlob.encode(.{ .id = 1, .active = true, .color = .red, .pos = .{ .x = 1, .y = 2 }, .tags = .{ 9, 9, 9 } }, &buf);
    var n: usize = 0;
    while (n < buf.len) : (n += 1) {
        try testing.expectError(error.TooShort, SampleBlob.decode(buf[0..n]));
    }
}

test "decode rejects wrong magic" {
    var buf: [SampleBlob.size]u8 = undefined;
    SampleBlob.encode(.{ .id = 1, .active = true, .color = .red, .pos = .{ .x = 1, .y = 2 }, .tags = .{ 9, 9, 9 } }, &buf);
    buf[0] = 'X';
    try testing.expectError(error.BadMagic, SampleBlob.decode(&buf));
}

test "decode rejects wrong version" {
    var buf: [SampleBlob.size]u8 = undefined;
    SampleBlob.encode(.{ .id = 1, .active = true, .color = .red, .pos = .{ .x = 1, .y = 2 }, .tags = .{ 9, 9, 9 } }, &buf);
    buf[4] = 2;
    buf[5] = 0;
    try testing.expectError(error.BadVersion, SampleBlob.decode(&buf));
}

test "decode rejects a corrupted checksum" {
    var buf: [SampleBlob.size]u8 = undefined;
    SampleBlob.encode(.{ .id = 1, .active = true, .color = .red, .pos = .{ .x = 1, .y = 2 }, .tags = .{ 9, 9, 9 } }, &buf);
    // Flip a payload byte without touching the trailer: checksum must catch it.
    buf[6] ^= 0xFF;
    try testing.expectError(error.BadChecksum, SampleBlob.decode(&buf));
}

test "fuzz decode never panics on arbitrary bytes" {
    try std.testing.fuzz({}, struct {
        fn f(_: void, s: *std.testing.Smith) !void {
            var buf: [512]u8 = undefined;
            const n = s.slice(&buf);
            _ = SampleBlob.decode(buf[0..n]) catch return;
        }
    }.f, .{});
}

test "mutated blobs are rejected, never accepted corrupt" {
    var good: [SampleBlob.size]u8 = undefined;
    SampleBlob.encode(std.mem.zeroes(Sample), &good);
    var prng: std.Random.DefaultPrng = .init(0xb10b);
    const r = prng.random();
    for (0..20_000) |_| {
        var m = good;
        const i = r.uintLessThan(usize, m.len);
        m[i] ^= r.intRangeAtMost(u8, 1, 255);
        const len = r.uintAtMost(usize, m.len);
        if (SampleBlob.decode(m[0..len])) |_| return error.TestUnexpectedResult else |_| {}
    }
}
