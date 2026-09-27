const rl = @import("raylib");

/// Move `camera.target` only once `target` leaves a dead-zone rectangle
/// centered in the viewport (as a fraction of screen size per axis in
/// `threshold`), then snap the zone's edge back onto it. Keeps the camera
/// still during small motion and catches up smoothly otherwise.
pub fn follow(camera: *rl.Camera2D, target: rl.Vector2, threshold: rl.Vector2) void {
    const display = rl.Vector2.init(
        @as(f32, @floatFromInt(rl.getRenderWidth())),
        @as(f32, @floatFromInt(rl.getRenderHeight())),
    );

    // Recenter the offset on the viewport every frame so a resize (or a
    // HiDPI zoom change) keeps the dead-zone centered rather than drifting.
    const half_min = rl.Vector2.one().subtract(threshold).scale(0.5).multiply(display);
    const half_max = rl.Vector2.one().add(threshold).scale(0.5).multiply(display);
    camera.offset = half_min;

    const world_min = rl.getScreenToWorld2D(half_min, camera.*);
    const world_max = rl.getScreenToWorld2D(half_max, camera.*);
    if (target.x < world_min.x) camera.target.x = target.x;
    if (target.y < world_min.y) camera.target.y = target.y;
    if (target.x > world_max.x) camera.target.x = world_min.x + (target.x - world_max.x);
    if (target.y > world_max.y) camera.target.y = world_min.y + (target.y - world_max.y);
}

/// Zoom factor that keeps one render pixel equal to one physical display
/// pixel on a HiDPI screen (`getRenderWidth`/`getRenderHeight` are the
/// framebuffer size; `getScreenWidth`/`getScreenHeight` are logical points).
pub fn hidpiZoom() f32 {
    const render_width: f32 = @floatFromInt(rl.getRenderWidth());
    const screen_width: f32 = @floatFromInt(rl.getScreenWidth());
    return render_width / screen_width;
}

// No tests: every function here calls into a live raylib context (a window
// and GL state), which doesn't exist under `zig test`. See RULES.md — rl/
// is compile-checked, not unit-tested.
