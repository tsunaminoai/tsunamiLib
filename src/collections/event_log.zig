/// EventLog(E) — an append / drain / clear-per-turn log for a poll + apply +
/// drain driver: the caller clears the log at the start of each `apply()` so
/// one turn's events never bleed into the next, then drains a slice (valid
/// until the next mutation) to narrate what happened. `clear` retains
/// capacity, so a whole run reuses one allocation.
const std = @import("std");

pub fn EventLog(comptime E: type) type {
    return struct {
        const Self = @This();

        events: std.ArrayList(E) = .empty,

        pub const empty: Self = .{ .events = .empty };

        pub fn deinit(self: *Self, gpa: std.mem.Allocator) void {
            self.events.deinit(gpa);
        }

        pub fn append(self: *Self, gpa: std.mem.Allocator, event: E) !void {
            try self.events.append(gpa, event);
        }

        /// The events accumulated since the last `clear()`. The slice is valid
        /// until the next mutation; frontends drain and narrate each turn.
        pub fn drain(self: *const Self) []const E {
            return self.events.items;
        }

        pub fn len(self: Self) usize {
            return self.events.items.len;
        }

        /// Drop all events, keeping capacity. Call at the start of each
        /// action so a turn's events don't leak into the next.
        pub fn clear(self: *Self) void {
            self.events.clearRetainingCapacity();
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const tst = std.testing;

const TestEvent = union(enum) {
    scored: struct { seat: u8, points: u8 },
    turn_passed: u8,
};

test "append then drain returns events in order" {
    var log: EventLog(TestEvent) = .empty;
    defer log.deinit(tst.allocator);

    try log.append(tst.allocator, .{ .scored = .{ .seat = 0, .points = 2 } });
    try log.append(tst.allocator, .{ .turn_passed = 1 });

    const drained = log.drain();
    try tst.expectEqual(@as(usize, 2), drained.len);
    try tst.expectEqual(@as(u8, 2), drained[0].scored.points);
    try tst.expectEqual(@as(u8, 1), drained[1].turn_passed);
}

test "clear empties the log but a fresh turn reuses it" {
    var log: EventLog(TestEvent) = .empty;
    defer log.deinit(tst.allocator);

    try log.append(tst.allocator, .{ .turn_passed = 0 });
    try tst.expectEqual(@as(usize, 1), log.len());

    // Start of the next turn: clear, then this turn's events accumulate fresh.
    log.clear();
    try tst.expectEqual(@as(usize, 0), log.len());
    try tst.expectEqual(@as(usize, 0), log.drain().len);

    try log.append(tst.allocator, .{ .scored = .{ .seat = 1, .points = 1 } });
    try tst.expectEqual(@as(usize, 1), log.drain().len);
}
