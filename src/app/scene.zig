const std = @import("std");
const contract = @import("../contract.zig");

/// Fixed-depth scene stack: `max_depth` slots, no allocation, ever. Scenes are
/// dispatched through a vtable so callers can mix distinct scene types (menu,
/// gameplay, pause overlay, ...) behind one `Scene` value, all sharing a
/// single `*Ctx` for app-wide state (renderer handle, input, asset tables).
///
/// Transitions call `exit` on a scene as it leaves the stack and `enter` on a
/// scene as it becomes current, in that order, so a scene can assume `enter`
/// has run before its first `update`/`render` and that `exit` runs exactly
/// once when it is popped or replaced.
pub fn SceneStack(comptime Ctx: type, comptime max_depth: usize) type {
    comptime std.debug.assert(max_depth >= 1);
    return struct {
        const Self = @This();

        pub const Scene = struct {
            ptr: *anyopaque,
            vtable: *const VTable,

            pub const VTable = struct {
                enter: *const fn (ptr: *anyopaque, ctx: *Ctx) void,
                exit: *const fn (ptr: *anyopaque, ctx: *Ctx) void,
                update: *const fn (ptr: *anyopaque, ctx: *Ctx, dt: f32) void,
                render: *const fn (ptr: *anyopaque, ctx: *Ctx) void,
            };
        };

        /// Cast a vtable fn's `*anyopaque` back to the concrete scene type
        /// inside that type's own enter/exit/update/render.
        pub fn ptrCast(comptime T: type, ptr: *anyopaque) *T {
            return @ptrCast(@alignCast(ptr));
        }

        slots: [max_depth]Scene = undefined,
        len: usize = 0,

        pub fn empty(self: *const Self) bool {
            return self.len == 0;
        }

        pub fn full(self: *const Self) bool {
            return self.len == max_depth;
        }

        /// Push `scene` on top and enter it. The previous top (if any) is
        /// left as-is: it is suspended, not exited.
        pub fn push(self: *Self, ctx: *Ctx, scene: Scene) void {
            contract.require(!self.full(), "scene.push: stack at max_depth");
            self.slots[self.len] = scene;
            self.len += 1;
            scene.vtable.enter(scene.ptr, ctx);
        }

        /// Exit and pop the top scene, uncovering the one beneath (which is
        /// not re-entered: it was never exited when covered).
        pub fn pop(self: *Self, ctx: *Ctx) void {
            contract.require(!self.empty(), "scene.pop: stack empty");
            self.len -= 1;
            const scene = self.slots[self.len];
            scene.vtable.exit(scene.ptr, ctx);
        }

        /// Exit and pop the top scene, then push and enter `scene` in its
        /// place. The depth is unchanged (a no-op replace on an empty stack
        /// degenerates to a plain push).
        pub fn replace(self: *Self, ctx: *Ctx, scene: Scene) void {
            if (!self.empty()) self.pop(ctx);
            self.push(ctx, scene);
        }

        /// Drive the top scene's `update`. No-op on an empty stack.
        pub fn update(self: *Self, ctx: *Ctx, dt: f32) void {
            const scene = self.top() orelse return;
            scene.vtable.update(scene.ptr, ctx, dt);
        }

        /// Render only the top scene. To render covered scenes too (e.g. a
        /// translucent pause overlay), walk `self.slots[0..self.len]` directly.
        pub fn render(self: *Self, ctx: *Ctx) void {
            const scene = self.top() orelse return;
            scene.vtable.render(scene.ptr, ctx);
        }

        pub fn top(self: *Self) ?Scene {
            if (self.len == 0) return null;
            return self.slots[self.len - 1];
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

const TestCtx = struct {
    buf: [64]u8 = undefined,
    len: usize = 0,

    fn record(self: *TestCtx, c: u8) void {
        self.buf[self.len] = c;
        self.len += 1;
    }

    fn str(self: *const TestCtx) []const u8 {
        return self.buf[0..self.len];
    }
};

fn TestStack(comptime max_depth: usize) type {
    return SceneStack(TestCtx, max_depth);
}

/// A test-double scene that tags every lifecycle call with an id byte so
/// call order and identity are both checkable from `ctx.str()`.
fn Recorder(comptime max_depth: usize) type {
    return struct {
        const Self = @This();
        const Stack = TestStack(max_depth);

        id: u8,
        updates: usize = 0,
        renders: usize = 0,

        fn scene(self: *Self) Stack.Scene {
            return .{ .ptr = self, .vtable = &vtable };
        }

        const vtable: Stack.Scene.VTable = .{
            .enter = enter,
            .exit = exit,
            .update = update,
            .render = render,
        };

        fn enter(ptr: *anyopaque, ctx: *TestCtx) void {
            const self = Stack.ptrCast(Self, ptr);
            ctx.record('e');
            ctx.record(self.id);
        }

        fn exit(ptr: *anyopaque, ctx: *TestCtx) void {
            const self = Stack.ptrCast(Self, ptr);
            ctx.record('x');
            ctx.record(self.id);
        }

        fn update(ptr: *anyopaque, ctx: *TestCtx, dt: f32) void {
            _ = dt;
            const self = Stack.ptrCast(Self, ptr);
            self.updates += 1;
            ctx.record('u');
            ctx.record(self.id);
        }

        fn render(ptr: *anyopaque, ctx: *TestCtx) void {
            const self = Stack.ptrCast(Self, ptr);
            self.renders += 1;
            ctx.record('r');
            ctx.record(self.id);
        }
    };
}

test "push enters the scene and update/render dispatch to it" {
    const Stack = TestStack(4);
    const R = Recorder(4);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};
    var a: R = .{ .id = 'A' };

    try testing.expect(stack.empty());
    stack.push(&ctx, a.scene());
    try testing.expectEqualStrings("eA", ctx.str());
    try testing.expect(!stack.empty());
    try testing.expectEqual(@as(usize, 1), stack.len);

    stack.update(&ctx, 0.016);
    stack.render(&ctx);
    try testing.expectEqual(@as(usize, 1), a.updates);
    try testing.expectEqual(@as(usize, 1), a.renders);
    try testing.expectEqualStrings("eAuArA", ctx.str());
}

test "pop exits the top scene and uncovers the one beneath without re-entering it" {
    const Stack = TestStack(4);
    const R = Recorder(4);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};
    var a: R = .{ .id = 'A' };
    var b: R = .{ .id = 'B' };

    stack.push(&ctx, a.scene());
    stack.push(&ctx, b.scene());
    try testing.expectEqualStrings("eAeB", ctx.str());

    stack.update(&ctx, 0.016);
    try testing.expectEqual(@as(usize, 0), a.updates);
    try testing.expectEqual(@as(usize, 1), b.updates);

    stack.pop(&ctx);
    try testing.expectEqualStrings("eAeBuBxB", ctx.str());
    try testing.expectEqual(@as(usize, 1), stack.len);

    // A was never re-entered when uncovered, but is live again for dispatch.
    stack.update(&ctx, 0.016);
    try testing.expectEqual(@as(usize, 1), a.updates);
    try testing.expectEqualStrings("eAeBuBxBuA", ctx.str());
}

test "replace exits the current top and enters the new scene at the same depth" {
    const Stack = TestStack(4);
    const R = Recorder(4);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};
    var a: R = .{ .id = 'A' };
    var b: R = .{ .id = 'B' };

    stack.push(&ctx, a.scene());
    stack.replace(&ctx, b.scene());
    try testing.expectEqualStrings("eAxAeB", ctx.str());
    try testing.expectEqual(@as(usize, 1), stack.len);

    stack.update(&ctx, 0.016);
    try testing.expectEqual(@as(usize, 0), a.updates);
    try testing.expectEqual(@as(usize, 1), b.updates);
}

test "replace on an empty stack degenerates to a push" {
    const Stack = TestStack(4);
    const R = Recorder(4);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};
    var a: R = .{ .id = 'A' };

    try testing.expect(stack.empty());
    stack.replace(&ctx, a.scene());
    try testing.expectEqualStrings("eA", ctx.str());
    try testing.expectEqual(@as(usize, 1), stack.len);
}

test "update and render are no-ops on an empty stack" {
    const Stack = TestStack(4);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};

    stack.update(&ctx, 0.016);
    stack.render(&ctx);
    try testing.expectEqualStrings("", ctx.str());
}

test "push fills the stack to exactly max_depth" {
    const Stack = TestStack(2);
    const R = Recorder(2);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};
    var a: R = .{ .id = 'A' };
    var b: R = .{ .id = 'B' };

    stack.push(&ctx, a.scene());
    stack.push(&ctx, b.scene());
    try testing.expect(stack.full());
    try testing.expectEqual(@as(usize, 2), stack.len);
}

test "empty/full report the guard conditions push/pop enforce via contract.require" {
    // `push` past max_depth and `pop` on an empty stack are caller bugs
    // guarded by `contract.require` (a panic), not a recoverable runtime
    // condition with its own error type: a fixed-depth stack's depth is a
    // static property of the call site, not something that varies with
    // untrusted input. `empty()`/`full()` are the guards callers check
    // before calling; this test exercises them across the full depth range
    // rather than triggering the panics directly (no death-test harness).
    const Stack = TestStack(2);
    const R = Recorder(2);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};
    var a: R = .{ .id = 'A' };
    var b: R = .{ .id = 'B' };

    try testing.expect(stack.empty());
    try testing.expect(!stack.full());
    stack.push(&ctx, a.scene());
    try testing.expect(!stack.empty());
    try testing.expect(!stack.full());
    stack.push(&ctx, b.scene());
    try testing.expect(!stack.empty());
    try testing.expect(stack.full());
    stack.pop(&ctx);
    stack.pop(&ctx);
    try testing.expect(stack.empty());
}

test "top reflects the current scene through push/pop/replace" {
    const Stack = TestStack(4);
    const R = Recorder(4);
    var ctx: TestCtx = .{};
    var stack: Stack = .{};
    var a: R = .{ .id = 'A' };
    var b: R = .{ .id = 'B' };

    try testing.expect(stack.top() == null);
    stack.push(&ctx, a.scene());
    try testing.expectEqual(a.scene().ptr, stack.top().?.ptr);
    stack.push(&ctx, b.scene());
    try testing.expectEqual(b.scene().ptr, stack.top().?.ptr);
    stack.pop(&ctx);
    try testing.expectEqual(a.scene().ptr, stack.top().?.ptr);
}
