//! Reed-Solomon erasure coding (k=4,m=2): encode, drop 2 of 6 segments,
//! reconstruct bit-exact. Also spreads a burst error across codewords
//! with a block interleaver so no single codeword sees more than its share.
const std = @import("std");
const ts = @import("tsunami");
const rs = ts.coding.rs;
const interleave = ts.coding.interleave;

const k = 4;
const m = 2;
const stripe = 32;

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;
    const gpa = init.gpa;

    var prng: std.Random.DefaultPrng = .init(0xE55);
    const rand = prng.random();

    var store: [k + m][stripe]u8 = undefined;
    for (0..k) |d| for (&store[d]) |*b| {
        b.* = rand.int(u8);
    };
    const orig = store;

    var coder = try rs.Coder.init(gpa, k, m);
    defer coder.deinit(gpa);

    var data: [k][]const u8 = undefined;
    for (0..k) |d| data[d] = &store[d];
    var parity: [m][]u8 = undefined;
    for (0..m) |p| parity[p] = &store[k + p];
    coder.encode(&data, &parity, stripe);

    // Drop 2 segments (one data, one parity) and clobber them to prove
    // reconstruction, not accidental survival.
    var present = [_]bool{true} ** (k + m);
    present[1] = false;
    present[k] = false;
    @memset(&store[1], 0xAA);
    @memset(&store[k], 0xAA);

    var segs: [k + m][]u8 = undefined;
    for (0..k + m) |j| segs[j] = &store[j];
    try coder.reconstruct(gpa, &segs, &present, stripe);

    var exact = true;
    for (0..k) |d| if (!std.mem.eql(u8, &orig[d], &store[d])) {
        exact = false;
    };

    try w.print("RS({d},{d}): dropped segments 1 (data) and {d} (parity)\n", .{ k + m, k, k });
    try w.print("reconstructed data bit-exact: {}\n", .{exact});

    // Interleaving: spread a contiguous burst of `depth` codeword-bits across
    // `width` codewords so each one absorbs exactly one corrupted bit.
    const depth = 8;
    const width = 64;
    const bits = try gpa.alloc(u1, depth * width);
    defer gpa.free(bits);
    @memset(bits, 0);
    const mixed = try gpa.alloc(u1, depth * width);
    defer gpa.free(mixed);
    interleave.interleave(u1, depth, width, bits, mixed);

    // Burst of `depth` consecutive errors in the transmitted (interleaved) stream.
    for (0..depth) |i| mixed[100 + i] = 1;

    const back = try gpa.alloc(u1, depth * width);
    defer gpa.free(back);
    interleave.deinterleave(u1, depth, width, mixed, back);

    var max_per_row: usize = 0;
    for (0..depth) |r| {
        var cnt: usize = 0;
        for (0..width) |c| cnt += back[r * width + c];
        max_per_row = @max(max_per_row, cnt);
    }
    try w.print("burst of {d} spread across {d} codewords: max {d} bit(s) in any one\n", .{ depth, depth, max_per_row });
    try w.flush();

    if (!exact or max_per_row != 1) return error.ExampleFailed;
}
