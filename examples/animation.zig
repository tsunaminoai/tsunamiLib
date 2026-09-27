//! One second at 60 fps, fixed dt: a Tween with two easings, a UTF-8
//! Typewriter reveal, a repeating Countdown, and a SceneStack push/pop.
const std = @import("std");
const ts = @import("tsunami");
const app = ts.app;

const dt: f32 = 1.0 / 60.0;
const frames = 60;

const Ctx = struct {
    log: [32]u8 = undefined,
    len: usize = 0,
    fn record(self: *Ctx, c: u8) void {
        self.log[self.len] = c;
        self.len += 1;
    }
};
const Stack = app.scene.SceneStack(Ctx, 4);

const Menu = struct {
    fn scene(self: *Menu) Stack.Scene {
        return .{ .ptr = self, .vtable = &vt };
    }
    const vt: Stack.Scene.VTable = .{ .enter = enter, .exit = exit, .update = noopDt, .render = noop };
    fn enter(_: *anyopaque, ctx: *Ctx) void {
        ctx.record('M');
    }
    fn exit(_: *anyopaque, ctx: *Ctx) void {
        ctx.record('m');
    }
    fn noop(_: *anyopaque, _: *Ctx) void {}
    fn noopDt(_: *anyopaque, _: *Ctx, _: f32) void {}
};

const Play = struct {
    fn scene(self: *Play) Stack.Scene {
        return .{ .ptr = self, .vtable = &vt };
    }
    const vt: Stack.Scene.VTable = .{ .enter = enter, .exit = exit, .update = Menu.noopDt, .render = Menu.noop };
    fn enter(_: *anyopaque, ctx: *Ctx) void {
        ctx.record('P');
    }
    fn exit(_: *anyopaque, ctx: *Ctx) void {
        ctx.record('p');
    }
};

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    var linear = app.tween.Tween(f32, app.tween.ease.linear).init(0, 100, 1.0);
    var back = app.tween.Tween(f32, app.tween.ease.backOut).init(0, 100, 1.0);
    var tw: app.typewriter.Typewriter = .{ .chars_per_second = 20 };
    const text = "hello, tsunami!";
    var beep: app.stopwatch.Countdown(f32) = .init(0.25, true);
    var beeps: usize = 0;

    var sw = app.stopwatch.Stopwatch.start(init.io);

    for (0..frames) |_| {
        linear.update(dt);
        back.update(dt);
        tw.update(text, dt);
        if (beep.update(dt)) beeps += 1;
    }
    const wall = sw.elapsed(init.io);

    var ctx: Ctx = .{};
    var stack: Stack = .{};
    var menu: Menu = .{};
    var play: Play = .{};
    stack.push(&ctx, menu.scene());
    stack.push(&ctx, play.scene());
    stack.pop(&ctx);
    stack.pop(&ctx);

    try w.print("linear @1s = {d:.0}, backOut @1s = {d:.1}\n", .{ linear.value(), back.value() });
    try w.print("typewriter: \"{s}\" (complete={})\n", .{ tw.revealed(text), tw.complete(text) });
    // f32 accumulation of 60 * (1/60) undershoots 1.0 slightly, so the 4th
    // 0.25s boundary lands just past the loop's end: 3 fires, not 4.
    try w.print("countdown fired {d} times in 1s (want 3)\n", .{beeps});
    try w.print("scene order: {s} (want MPpm)\n", .{ctx.log[0..ctx.len]});
    try w.print("60-frame loop took {d}ms wall time\n", .{wall.toMilliseconds()});

    const ok = linear.value() > 99.9 and back.value() > 99.9 and tw.complete(text) and beeps == 3 and
        std.mem.eql(u8, ctx.log[0..ctx.len], "MPpm");
    try w.flush();
    if (!ok) return error.ExampleFailed;
}
