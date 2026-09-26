/// Id(tag) and IdSlab: type-safe distinct IDs and their backing store. IDs
/// with different tags are distinct types, so you can't mix up IDs from
/// different stores without an explicit cast.
const std = @import("std");

pub fn Id(comptime tag: []const u8) type {
    return struct {
        const Self = @This();

        value: u32,

        pub const invalid: Self = .{ .value = std.math.maxInt(u32) };

        pub fn from(v: u32) Self {
            return .{ .value = v };
        }

        pub fn index(self: Self) u32 {
            return self.value;
        }

        pub fn format(self: Self, comptime fmt: []const u8, _: std.fmt.FormatOptions, w: *std.Io.Writer) std.Io.Writer.Error!void {
            _ = fmt;
            try w.print("{s}({})", .{ tag, self.value });
        }
    };
}

pub fn IdSlab(comptime T: type, comptime tag: []const u8) type {
    return struct {
        const Self = @This();
        pub const IdType = Id(tag);

        items: std.ArrayList(T) = .empty,

        pub const empty: Self = .{ .items = .empty };

        pub fn deinit(self: *Self, alloc: std.mem.Allocator) void {
            self.items.deinit(alloc);
        }

        pub fn insert(self: *Self, alloc: std.mem.Allocator, value: T) !IdType {
            const idx = self.items.items.len;
            try self.items.append(alloc, value);
            return IdType.from(@intCast(idx));
        }

        pub fn get(self: Self, id: IdType) ?T {
            const idx = id.index();
            if (idx >= self.items.items.len) return null;
            return self.items.items[idx];
        }

        pub fn getPtr(self: *Self, id: IdType) ?*T {
            const idx = id.index();
            if (idx >= self.items.items.len) return null;
            return &self.items.items[idx];
        }

        pub fn set(self: *Self, id: IdType, value: T) !void {
            const ptr = self.getPtr(id) orelse return error.InvalidId;
            ptr.* = value;
        }

        pub fn len(self: Self) usize {
            return self.items.items.len;
        }

        /// Clears all items but retains capacity.
        pub fn clear(self: *Self) void {
            self.items.clearRetainingCapacity();
        }

        pub fn clone(self: Self, alloc: std.mem.Allocator) !Self {
            return .{ .items = try self.items.clone(alloc) };
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const tst = std.testing;

test "Id types are distinct by tag" {
    const UserId = Id("user");
    const RoleId = Id("role");

    const user1 = UserId.from(42);
    const role1 = RoleId.from(42);

    // Both have the same underlying value but different types.
    try tst.expectEqual(user1.index(), role1.index());
    // The types are different (compile-time check in real code).
    try tst.expectEqual(@TypeOf(user1) == @TypeOf(role1), false);
}

test "IdSlab insert and get" {
    var slab: IdSlab(u32, "test") = .empty;
    defer slab.deinit(tst.allocator);

    const id1 = try slab.insert(tst.allocator, 100);
    const id2 = try slab.insert(tst.allocator, 200);

    try tst.expectEqual(@as(u32, 100), slab.get(id1).?);
    try tst.expectEqual(@as(u32, 200), slab.get(id2).?);

    // Out of range returns null.
    try tst.expectEqual(@as(?u32, null), slab.get(.{ .value = 999 }));
}

test "IdSlab set and getPtr" {
    var slab: IdSlab(i32, "test") = .empty;
    defer slab.deinit(tst.allocator);

    const id = try slab.insert(tst.allocator, 10);
    try slab.set(id, 20);
    try tst.expectEqual(@as(i32, 20), slab.get(id).?);

    // Modify via getPtr.
    if (slab.getPtr(id)) |ptr| {
        ptr.* += 15;
    }
    try tst.expectEqual(@as(i32, 35), slab.get(id).?);
}

test "IdSlab len, clear, and clone" {
    var slab: IdSlab(u8, "test") = .empty;
    defer slab.deinit(tst.allocator);

    _ = try slab.insert(tst.allocator, 1);
    _ = try slab.insert(tst.allocator, 2);
    try tst.expectEqual(@as(usize, 2), slab.len());

    slab.clear();
    try tst.expectEqual(@as(usize, 0), slab.len());

    // Clone an empty slab.
    var cloned = try slab.clone(tst.allocator);
    defer cloned.deinit(tst.allocator);
    try tst.expectEqual(@as(usize, 0), cloned.len());
}

test "IdSlab clone copies data" {
    var original: IdSlab(u32, "test") = .empty;
    defer original.deinit(tst.allocator);

    const id1 = try original.insert(tst.allocator, 100);
    const id2 = try original.insert(tst.allocator, 200);

    var cloned = try original.clone(tst.allocator);
    defer cloned.deinit(tst.allocator);

    try tst.expectEqual(@as(u32, 100), cloned.get(id1).?);
    try tst.expectEqual(@as(u32, 200), cloned.get(id2).?);

    // Modifying clone doesn't affect original.
    try cloned.set(id1, 999);
    try tst.expectEqual(@as(u32, 100), original.get(id1).?);
    try tst.expectEqual(@as(u32, 999), cloned.get(id1).?);
}
