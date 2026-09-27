/// Pile(T) — an ordered container of cards (or any value type) for card games.
/// `append` never sorts implicitly (order is data — pegging scores depend on
/// play order, not rank order); ask for order explicitly via `insertSorted`
/// or `sort`. Unmanaged: stores no allocator, so it can be embedded by value.
const std = @import("std");

pub fn Pile(comptime T: type) type {
    return struct {
        const Self = @This();

        list: std.ArrayList(T) = .empty,

        /// An empty pile with no backing allocation. Embed as
        /// `pile: Pile(Card) = .empty` or assign at init.
        pub const empty: Self = .{ .list = .empty };

        pub fn deinit(self: *Self, gpa: std.mem.Allocator) void {
            self.list.deinit(gpa);
        }

        // ── size & access ────────────────────────────────────────────────────

        pub fn len(self: Self) usize {
            return self.list.items.len;
        }

        /// Mutable view of the contents in order. Handy for scoring passes that
        /// iterate the whole pile; do not retain across a mutation.
        pub fn items(self: Self) []T {
            return self.list.items;
        }

        /// Bounds-safe read: returns null for `idx >= len`.
        pub fn get(self: Self, idx: usize) ?T {
            if (idx >= self.len()) return null;
            return self.list.items[idx];
        }

        /// The last element without removing it, or null if empty.
        pub fn peek(self: Self) ?T {
            if (self.len() == 0) return null;
            return self.list.items[self.len() - 1];
        }

        // ── insertion ────────────────────────────────────────────────────────

        /// Append preserving insertion order. NO implicit sort — use
        /// `insertSorted` or `sort` if order is wanted.
        pub fn append(self: *Self, gpa: std.mem.Allocator, item: T) !void {
            try self.list.append(gpa, item);
        }

        /// Append several items in order (no sort).
        pub fn appendSlice(self: *Self, gpa: std.mem.Allocator, slice: []const T) !void {
            try self.list.appendSlice(gpa, slice);
        }

        /// Insert keeping the pile ordered by `lessThan`. `context`/`lessThan`
        /// mirror the `std.mem.sort` convention. Linear insertion; card-game
        /// piles are tiny.
        pub fn insertSorted(
            self: *Self,
            gpa: std.mem.Allocator,
            item: T,
            context: anytype,
            comptime lessThan: fn (@TypeOf(context), T, T) bool,
        ) !void {
            try self.list.ensureUnusedCapacity(gpa, 1);
            var i: usize = self.len();
            while (i > 0 and lessThan(context, item, self.list.items[i - 1])) : (i -= 1) {}
            self.list.insertAssumeCapacity(i, item);
        }

        /// Sort in place by `lessThan`. Same convention as `std.mem.sort`.
        pub fn sort(
            self: *Self,
            context: anytype,
            comptime lessThan: fn (@TypeOf(context), T, T) bool,
        ) void {
            std.mem.sort(T, self.list.items, context, lessThan);
        }

        // ── removal ──────────────────────────────────────────────────────────

        /// Remove and return the top element (the last appended) — a deck draw.
        /// Null if empty.
        pub fn draw(self: *Self) ?T {
            return self.list.pop();
        }

        /// Bounds-safe ordered remove: null for `idx >= len`, else removes at
        /// `idx` (shifting the tail down) and returns the element.
        pub fn remove(self: *Self, idx: usize) ?T {
            if (idx >= self.len()) return null;
            return self.list.orderedRemove(idx);
        }

        // ── shuffle ──────────────────────────────────────────────────────────

        /// Shuffle in place. `rand` is passed per call, never stored — see
        /// rand/rng.zig for why a stored interface pointer would dangle.
        pub fn shuffle(self: Self, rand: std.Random) void {
            rand.shuffle(T, self.list.items);
        }

        // ── clearing ─────────────────────────────────────────────────────────

        /// Empty the pile but keep the backing capacity for reuse between
        /// rounds. Refilling without clearing first silently accumulates.
        pub fn clear(self: *Self) void {
            self.list.clearRetainingCapacity();
        }

        /// Empty the pile and release the backing capacity.
        pub fn clearAndFree(self: *Self, gpa: std.mem.Allocator) void {
            self.list.clearAndFree(gpa);
        }

        /// A deep copy. Caller owns the result and must `deinit` it.
        pub fn clone(self: Self, gpa: std.mem.Allocator) !Self {
            return .{ .list = try self.list.clone(gpa) };
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const tst = std.testing;

const TestCard = struct {
    rank: u8,
    fn lessThan(_: void, a: TestCard, b: TestCard) bool {
        return a.rank < b.rank;
    }
};

test "append preserves insertion order (no implicit sort)" {
    var p: Pile(TestCard) = .empty;
    defer p.deinit(tst.allocator);

    try p.append(tst.allocator, .{ .rank = 5 });
    try p.append(tst.allocator, .{ .rank = 1 });
    try p.append(tst.allocator, .{ .rank = 5 });

    // Exactly the order played — a sort-on-insert pile would have put the two
    // 5s next to each other and mis-scored a pegging pair.
    try tst.expectEqual(@as(usize, 3), p.len());
    try tst.expectEqual(@as(u8, 5), p.get(0).?.rank);
    try tst.expectEqual(@as(u8, 1), p.get(1).?.rank);
    try tst.expectEqual(@as(u8, 5), p.get(2).?.rank);
}

test "get and remove are bounds-safe" {
    var p: Pile(TestCard) = .empty;
    defer p.deinit(tst.allocator);

    try p.append(tst.allocator, .{ .rank = 1 });
    try p.append(tst.allocator, .{ .rank = 2 });

    // idx == len must be null (the classic off-by-one that read one past end).
    try tst.expect(p.get(p.len()) == null);
    try tst.expect(p.get(999) == null);
    try tst.expect(p.remove(p.len()) == null);

    try tst.expect(p.get(0) != null);
    try tst.expect(p.get(p.len() - 1) != null);

    const removed = p.remove(0).?;
    try tst.expectEqual(@as(u8, 1), removed.rank);
    try tst.expectEqual(@as(usize, 1), p.len());
    try tst.expectEqual(@as(u8, 2), p.get(0).?.rank);
}

test Pile {
    var p: Pile(TestCard) = .empty;
    defer p.deinit(tst.allocator);

    try tst.expect(p.peek() == null);
    try tst.expect(p.draw() == null);

    try p.append(tst.allocator, .{ .rank = 7 });
    try p.append(tst.allocator, .{ .rank = 9 });

    try tst.expectEqual(@as(u8, 9), p.peek().?.rank); // does not remove
    try tst.expectEqual(@as(usize, 2), p.len());
    try tst.expectEqual(@as(u8, 9), p.draw().?.rank); // removes the top
    try tst.expectEqual(@as(usize, 1), p.len());
}

test "insertSorted and sort produce the same order for unique keys" {
    var a: Pile(TestCard) = .empty;
    defer a.deinit(tst.allocator);
    var b: Pile(TestCard) = .empty;
    defer b.deinit(tst.allocator);

    const feed = [_]u8{ 3, 1, 4, 2, 5 };
    for (feed) |r| {
        try a.insertSorted(tst.allocator, .{ .rank = r }, {}, TestCard.lessThan);
        try b.append(tst.allocator, .{ .rank = r });
    }
    b.sort({}, TestCard.lessThan);

    try tst.expectEqual(a.len(), b.len());
    for (0..a.len()) |i| {
        try tst.expectEqual(a.get(i).?.rank, b.get(i).?.rank);
        try tst.expectEqual(@as(u8, @intCast(i + 1)), a.get(i).?.rank);
    }
}

test "shuffle is a permutation and deterministic for a seed" {
    var p: Pile(u8) = .empty;
    defer p.deinit(tst.allocator);
    for (0..52) |i| try p.append(tst.allocator, @intCast(i));

    var r1 = std.Random.DefaultPrng.init(0xC0FFEE);
    p.shuffle(r1.random());

    // Still a permutation of 0..52 (no loss, no duplication — the reuse/clear
    // discipline plus shuffle must not change the multiset).
    var seen = [_]bool{false} ** 52;
    for (p.items()) |v| {
        try tst.expect(!seen[v]);
        seen[v] = true;
    }

    // Same seed reproduces the same order.
    var q: Pile(u8) = .empty;
    defer q.deinit(tst.allocator);
    for (0..52) |i| try q.append(tst.allocator, @intCast(i));
    var r2 = std.Random.DefaultPrng.init(0xC0FFEE);
    q.shuffle(r2.random());
    try tst.expectEqualSlices(u8, p.items(), q.items());
}

test "clear keeps capacity; refill does not accumulate" {
    var p: Pile(u8) = .empty;
    defer p.deinit(tst.allocator);

    for (0..52) |i| try p.append(tst.allocator, @intCast(i));
    for (0..13) |_| _ = p.draw();
    try tst.expectEqual(@as(usize, 39), p.len());

    // Reuse: clear, then refill to exactly 52 (not 39 + 52 = 91).
    p.clear();
    try tst.expectEqual(@as(usize, 0), p.len());
    for (0..52) |i| try p.append(tst.allocator, @intCast(i));
    try tst.expectEqual(@as(usize, 52), p.len());
}

test "clone is independent" {
    var p: Pile(u8) = .empty;
    defer p.deinit(tst.allocator);
    try p.append(tst.allocator, 1);
    try p.append(tst.allocator, 2);

    var c = try p.clone(tst.allocator);
    defer c.deinit(tst.allocator);

    _ = p.draw();
    try tst.expectEqual(@as(usize, 1), p.len());
    try tst.expectEqual(@as(usize, 2), c.len()); // clone unaffected
}
