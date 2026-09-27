const std = @import("std");
const rl = @import("raylib");
const ts = @import("tsunami");

/// Draw the currently-revealed prefix of `text` (see `ts.app.typewriter`),
/// wrapping is the caller's problem — this draws one line via raylib's
/// default font. `tw` must have been `update`d against this same `text`.
pub fn typewriter(tw: *const ts.app.typewriter.Typewriter, text: [:0]const u8, x: i32, y: i32, size: i32, color: rl.Color) void {
    const revealed = tw.revealed(text);
    // revealed() is a slice of `text`; re-terminate so drawText gets a
    // sentinel without a copy when it already ends the string, or with one
    // when it's a strict prefix.
    if (revealed.len == text.len) {
        rl.drawText(text, x, y, size, color);
    } else {
        var buf: [512]u8 = undefined;
        const n = @min(revealed.len, buf.len - 1);
        @memcpy(buf[0..n], revealed[0..n]);
        buf[n] = 0;
        rl.drawText(buf[0..n :0], x, y, size, color);
    }
}

/// Draw a vertical stack of buttons and report which one (if any) was
/// clicked this frame. `labels[i]` is drawn at `rectFn(i)`.
pub fn menu(labels: []const [:0]const u8, rectFn: *const fn (usize) rl.Rectangle, mouse: rl.Vector2, clicked: bool, drawButton: *const fn (rl.Rectangle, [:0]const u8, bool, bool) void) ?usize {
    var picked: ?usize = null;
    for (labels, 0..) |label, i| {
        const rect = rectFn(i);
        const hovered = rl.checkCollisionPointRec(mouse, rect);
        drawButton(rect, label, hovered, true);
        if (hovered and clicked) picked = i;
    }
    return picked;
}

/// Draw a cubic Bezier wire between two screen-space points, control points
/// pulled horizontally outward from each endpoint (the standard "wire" look
/// for a node-graph editor). Uses `ts.math.geom.bezierPoint` for the curve
/// math so the geometry is shared with hit-testing.
pub fn wire(from: rl.Vector2, to: rl.Vector2, thickness: f32, color: rl.Color) void {
    const V = @Vector(2, f32);
    const spread = @max(@abs(to.x - from.x) * 0.4 + 30, 40.0);
    const p0: V = .{ from.x, from.y };
    const p1: V = .{ from.x + spread, from.y };
    const p2: V = .{ to.x - spread, to.y };
    const p3: V = .{ to.x, to.y };

    const steps: usize = 24;
    var prev = from;
    for (1..steps + 1) |si| {
        const t: f32 = @as(f32, @floatFromInt(si)) / @as(f32, @floatFromInt(steps));
        const p = ts.math.geom.bezierPoint(V, t, p0, p1, p2, p3);
        const cur: rl.Vector2 = .{ .x = p[0], .y = p[1] };
        rl.drawLineEx(prev, cur, thickness, color);
        prev = cur;
    }
}

/// Whether `mouse` is within `thresh` px of the wire drawn by `wire` between
/// `from` and `to`. Shares the same control-point placement so a hit test
/// always matches what got drawn.
pub fn wireHitTest(from: rl.Vector2, to: rl.Vector2, mouse: rl.Vector2, thresh: f32) bool {
    const V = @Vector(2, f32);
    const spread = @max(@abs(to.x - from.x) * 0.4 + 30, 40.0);
    const p0: V = .{ from.x, from.y };
    const p1: V = .{ from.x + spread, from.y };
    const p2: V = .{ to.x - spread, to.y };
    const p3: V = .{ to.x, to.y };
    const m: V = .{ mouse.x, mouse.y };
    return ts.math.geom.distPointBezier(V, 12, m, p0, p1, p2, p3) <= thresh;
}

// No tests: every function here draws through raylib (drawText/drawLineEx/
// checkCollisionPointRec need a live GL context). See RULES.md — rl/ is
// compile-checked, not unit-tested. The bezier math itself is tested in
// ts.math.geom.
