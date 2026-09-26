const rl = @import("raylib");

/// Semantic screen-region -> pixel-rectangle mapping, generic over a
/// caller-supplied comptime `Region` enum. `rectFn` is `fn(Region) Rectangle`
/// — normally a big `switch` living beside the app's own layout constants —
/// so this module owns only the lookup contract, not any app's geometry.
pub fn Layout(comptime Region: type) type {
    return struct {
        rectFn: *const fn (Region) rl.Rectangle,

        const Self = @This();

        pub fn init(rectFn: *const fn (Region) rl.Rectangle) Self {
            return .{ .rectFn = rectFn };
        }

        pub fn rect(self: Self, region: Region) rl.Rectangle {
            return self.rectFn(region);
        }

        /// Whether `point` falls inside `region`'s rectangle.
        pub fn contains(self: Self, region: Region, point: rl.Vector2) bool {
            return rl.checkCollisionPointRec(point, self.rect(region));
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────
// Only the pure lookup/dispatch logic is testable without a raylib context
// (checkCollisionPointRec is a trivial AABB test the real header implements
// without touching GL state, but we keep this module's own tests to the
// parts that don't need a window at all — the rect dispatch).

const testing = @import("std").testing;

const TestRegion = enum { header, body, footer };

fn testRect(r: TestRegion) rl.Rectangle {
    return switch (r) {
        .header => .{ .x = 0, .y = 0, .width = 100, .height = 20 },
        .body => .{ .x = 0, .y = 20, .width = 100, .height = 60 },
        .footer => .{ .x = 0, .y = 80, .width = 100, .height = 20 },
    };
}

test "Layout dispatches region to its rectangle" {
    const L = Layout(TestRegion).init(&testRect);
    const body = L.rect(.body);
    try testing.expectEqual(@as(f32, 20), body.y);
    try testing.expectEqual(@as(f32, 60), body.height);
}
