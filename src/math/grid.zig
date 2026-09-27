const std = @import("std");
const contract = @import("../contract.zig");

/// 2D integer vector: cell coordinates, offsets, directions.
pub fn Pos(comptime I: type) type {
    comptime std.debug.assert(@typeInfo(I) == .int);
    return struct {
        x: I = 0,
        y: I = 0,

        const Self = @This();

        pub fn init(x: I, y: I) Self {
            return .{ .x = x, .y = y };
        }

        pub fn add(a: Self, b: Self) Self {
            return .{ .x = a.x + b.x, .y = a.y + b.y };
        }

        pub fn sub(a: Self, b: Self) Self {
            return .{ .x = a.x - b.x, .y = a.y - b.y };
        }

        pub fn scale(a: Self, s: I) Self {
            return .{ .x = a.x * s, .y = a.y * s };
        }

        pub fn manhattan(a: Self, b: Self) I {
            return @intCast(@abs(a.x - b.x) + @abs(a.y - b.y));
        }

        /// Rotate 90° clockwise in a y-down (screen) frame: (x, y) -> (-y, x).
        pub fn rotate90(a: Self) Self {
            return .{ .x = -a.y, .y = a.x };
        }
    };
}

/// Dense row-major 2D grid over a caller-owned or gpa-allocated slice.
/// `index = y*width + x`.
pub fn Grid(comptime T: type) type {
    return struct {
        data: []T,
        width: usize,
        height: usize,
        owned: bool = false,

        const Self = @This();
        pub const P = Pos(i32);

        /// Wrap a caller-owned buffer; `buf.len` must equal `width*height`.
        pub fn initBuffer(buf: []T, width: usize, height: usize) Self {
            contract.require(buf.len == width * height, "Grid.initBuffer: buf.len mismatch");
            return .{ .data = buf, .width = width, .height = height };
        }

        /// Allocate a `width*height` buffer, all cells set to `fill`.
        pub fn initAlloc(gpa: std.mem.Allocator, width: usize, height: usize, fill: T) !Self {
            const buf = try gpa.alloc(T, width * height);
            @memset(buf, fill);
            return .{ .data = buf, .width = width, .height = height, .owned = true };
        }

        pub fn deinit(self: *Self, gpa: std.mem.Allocator) void {
            if (self.owned) gpa.free(self.data);
            self.* = undefined;
        }

        pub inline fn inBounds(self: Self, x: i32, y: i32) bool {
            return x >= 0 and y >= 0 and x < @as(i32, @intCast(self.width)) and y < @as(i32, @intCast(self.height));
        }

        inline fn index(self: Self, x: usize, y: usize) usize {
            return y * self.width + x;
        }

        pub fn get(self: Self, x: i32, y: i32) ?T {
            if (!self.inBounds(x, y)) return null;
            return self.data[self.index(@intCast(x), @intCast(y))];
        }

        pub fn set(self: *Self, x: i32, y: i32, v: T) void {
            contract.require(self.inBounds(x, y), "Grid.set: out of bounds");
            self.data[self.index(@intCast(x), @intCast(y))] = v;
        }

        /// Slice of row `y`, or null if out of range.
        pub fn row(self: Self, y: i32) ?[]T {
            if (y < 0 or y >= @as(i32, @intCast(self.height))) return null;
            const base = self.index(0, @intCast(y));
            return self.data[base .. base + self.width];
        }

        pub const Neighbor = struct { pos: P, value: T };

        const offsets4 = [_]P{ .{ .x = 0, .y = -1 }, .{ .x = 1, .y = 0 }, .{ .x = 0, .y = 1 }, .{ .x = -1, .y = 0 } };
        const offsets8 = [_]P{
            .{ .x = 0, .y = -1 }, .{ .x = 1, .y = -1 }, .{ .x = 1, .y = 0 },  .{ .x = 1, .y = 1 },
            .{ .x = 0, .y = 1 },  .{ .x = -1, .y = 1 }, .{ .x = -1, .y = 0 }, .{ .x = -1, .y = -1 },
        };

        fn NeighborIter(comptime offsets: []const P) type {
            return struct {
                grid: *const Self,
                center: P,
                i: usize = 0,

                pub fn next(it: *@This()) ?Neighbor {
                    while (it.i < offsets.len) {
                        const off = offsets[it.i];
                        it.i += 1;
                        const p = it.center.add(off);
                        if (it.grid.get(p.x, p.y)) |v| return .{ .pos = p, .value = v };
                    }
                    return null;
                }
            };
        }

        pub fn neighbors4(self: *const Self, x: i32, y: i32) NeighborIter(&offsets4) {
            return .{ .grid = self, .center = .{ .x = x, .y = y } };
        }

        pub fn neighbors8(self: *const Self, x: i32, y: i32) NeighborIter(&offsets8) {
            return .{ .grid = self, .center = .{ .x = x, .y = y } };
        }
    };
}

/// Dense 3D grid extents, Fortran memory order (i fastest), 0-based:
/// `index = i + ni*(j + nj*k)`. Pure index arithmetic — no allocation, no data.
pub const Grid3 = struct {
    ni: usize,
    nj: usize,
    nk: usize,

    pub inline fn at(g: Grid3, i: usize, j: usize, k: usize) usize {
        return i + g.ni * (j + g.nj * k);
    }

    pub inline fn len(g: Grid3) usize {
        return g.ni * g.nj * g.nk;
    }

    pub inline fn column(g: Grid3, i: usize, j: usize) Column {
        return .{ .base = i + g.ni * j, .stride = g.ni * g.nj };
    }
};

/// Strided 1D view into a flat slice: element k is at `base + k*stride`.
pub const Column = struct {
    base: usize,
    stride: usize,

    pub inline fn at(c: Column, k: usize) usize {
        return c.base + k * c.stride;
    }
};

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "Pos add/sub/scale/manhattan/rotate90" {
    const P = Pos(i32);
    const a = P.init(2, 3);
    const b = P.init(-1, 5);
    try testing.expectEqual(P.init(1, 8), a.add(b));
    try testing.expectEqual(P.init(3, -2), a.sub(b));
    try testing.expectEqual(P.init(4, 6), a.scale(2));
    try testing.expectEqual(@as(i32, 5), a.manhattan(b));
    try testing.expectEqual(P.init(-3, 2), P.init(2, 3).rotate90());
}

test Grid {
    var buf: [12]u8 = undefined;
    var g = Grid(u8).initBuffer(&buf, 4, 3);
    try testing.expect(g.inBounds(0, 0));
    try testing.expect(!g.inBounds(4, 0));
    try testing.expect(!g.inBounds(-1, 0));
    g.set(2, 1, 42);
    try testing.expectEqual(@as(?u8, 42), g.get(2, 1));
    try testing.expectEqual(@as(?u8, null), g.get(4, 4));
}

test "Grid initAlloc/deinit and row slicing" {
    var g = try Grid(i32).initAlloc(testing.allocator, 5, 2, -1);
    defer g.deinit(testing.allocator);
    for (g.row(0).?) |v| try testing.expectEqual(@as(i32, -1), v);
    try testing.expect(g.row(2) == null);
    g.set(3, 1, 7);
    try testing.expectEqual(@as(i32, 7), g.row(1).?[3]);
}

test "neighbors4 skips out-of-bounds" {
    var buf: [9]u8 = .{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    const g = Grid(u8).initBuffer(&buf, 3, 3);
    var it = g.neighbors4(0, 0);
    var count: usize = 0;
    var sum: usize = 0;
    while (it.next()) |n| {
        count += 1;
        sum += n.value;
    }
    try testing.expectEqual(@as(usize, 2), count); // only right (1) and down (3)
    try testing.expectEqual(@as(usize, 1 + 3), sum);
}

test "neighbors8 from center sees all 8" {
    var buf: [9]u8 = .{ 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    const g = Grid(u8).initBuffer(&buf, 3, 3);
    var it = g.neighbors8(1, 1);
    var count: usize = 0;
    while (it.next()) |_| count += 1;
    try testing.expectEqual(@as(usize, 8), count);
}

test "Grid3 indexing is Fortran order, i fastest" {
    const g: Grid3 = .{ .ni = 3, .nj = 4, .nk = 5 };
    try testing.expectEqual(@as(usize, 60), g.len());
    try testing.expectEqual(@as(usize, 0), g.at(0, 0, 0));
    try testing.expectEqual(@as(usize, 1), g.at(1, 0, 0));
    try testing.expectEqual(@as(usize, 3), g.at(0, 1, 0));
    try testing.expectEqual(@as(usize, 12), g.at(0, 0, 1));
    try testing.expectEqual(g.len() - 1, g.at(2, 3, 4));
}

test "Grid3.column addresses the same cells as at()" {
    const g: Grid3 = .{ .ni = 7, .nj = 3, .nk = 6 };
    const col = g.column(4, 2);
    for (0..g.nk) |k| {
        try testing.expectEqual(g.at(4, 2, k), col.at(k));
    }
}
