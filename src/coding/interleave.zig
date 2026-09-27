//! Block interleaver: the transpose spreads a burst of B stream symbols
//! across `depth` codewords as exactly B/depth each (B a multiple of depth).
const std = @import("std");
const contract = @import("../contract.zig");

/// out[c·depth + r] = in[r·width + c]; both slices are depth·width long.
pub fn interleave(comptime T: type, depth: usize, width: usize, in: []const T, out: []T) void {
    contract.require(in.len == depth * width and out.len == in.len, "interleave: len != depth*width");
    for (0..depth) |r| {
        for (0..width) |c| {
            out[c * depth + r] = in[r * width + c];
        }
    }
}

/// Inverse of `interleave`.
pub fn deinterleave(comptime T: type, depth: usize, width: usize, in: []const T, out: []T) void {
    contract.require(in.len == depth * width and out.len == in.len, "interleave: len != depth*width");
    for (0..depth) |r| {
        for (0..width) |c| {
            out[r * width + c] = in[c * depth + r];
        }
    }
}

// ── Tests ────────────────────────────────────────────────────────────────

test interleave {
    const alloc = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(21);
    const rand = prng.random();
    const depth = 16;
    const width = 1944;

    const bits = try alloc.alloc(u1, depth * width);
    defer alloc.free(bits);
    for (bits) |*b| b.* = @intCast(rand.int(u1));
    const mixed = try alloc.alloc(u1, depth * width);
    defer alloc.free(mixed);
    const back = try alloc.alloc(u1, depth * width);
    defer alloc.free(back);
    interleave(u1, depth, width, bits, mixed);
    deinterleave(u1, depth, width, mixed, back);
    try std.testing.expectEqualSlices(u1, bits, back);

    const llrs = try alloc.alloc(f32, depth * width);
    defer alloc.free(llrs);
    for (llrs) |*v| v.* = rand.floatNorm(f32);
    const mixed_f = try alloc.alloc(f32, depth * width);
    defer alloc.free(mixed_f);
    const back_f = try alloc.alloc(f32, depth * width);
    defer alloc.free(back_f);
    interleave(f32, depth, width, llrs, mixed_f);
    deinterleave(f32, depth, width, mixed_f, back_f);
    try std.testing.expectEqualSlices(f32, llrs, back_f);
}

test "576-bit burst lands exactly 36 bits in every codeword" {
    const alloc = std.testing.allocator;
    const depth = 16;
    const width = 1944;
    const mixed = try alloc.alloc(u1, depth * width);
    defer alloc.free(mixed);
    @memset(mixed, 0);
    // burst at an arbitrary unaligned position in the INTERLEAVED stream
    const burst_start = 12345;
    for (burst_start..burst_start + 576) |i| mixed[i] = 1;
    const back = try alloc.alloc(u1, depth * width);
    defer alloc.free(back);
    deinterleave(u1, depth, width, mixed, back);
    for (0..depth) |r| {
        var cnt: usize = 0;
        for (0..width) |c| cnt += back[r * width + c];
        try std.testing.expectEqual(@as(usize, 36), cnt);
    }
}
