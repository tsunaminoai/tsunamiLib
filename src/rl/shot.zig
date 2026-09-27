const std = @import("std");
const rl = @import("raylib");

/// One named frame to capture: `draw` runs with the render texture already
/// bound (inside `beginTextureMode`/`endTextureMode`), and its output is
/// exported to `<out_dir>/<name>.png`.
pub const Shot = struct {
    name: []const u8,
    draw: *const fn () void,
};

/// Options gathering everything env-var driven, so a caller can either read
/// them from the process environment (`fromEnviron`) or build them directly
/// for a test.
pub const Options = struct {
    /// Whether the harness should run at all (e.g. gated on `APP_SHOT=1`).
    enabled: bool,
    /// Output directory, relative to CWD (raylib's `exportImage` writes
    /// relative paths); the caller must ensure it exists.
    out_dir: []const u8 = "runs/shots",
};

/// Build `Options` from an environ map: `enabled` from `enabled_var`
/// (any non-empty value counts as enabled), `out_dir` from `dir_var` when
/// set. Takes the map rather than calling getenv so it's testable without a
/// real process environment.
pub fn optionsFromEnviron(map: *const std.process.Environ.Map, enabled_var: []const u8, dir_var: []const u8) Options {
    var opts: Options = .{ .enabled = false };
    if (map.get(enabled_var)) |v| opts.enabled = v.len > 0;
    if (map.get(dir_var)) |d| opts.out_dir = d;
    return opts;
}

/// Render every shot in `shots` into an off-screen `width`x`height` render
/// texture and export each as a PNG under `opts.out_dir`. Uses `exportImage`
/// off the render texture, NOT `takeScreenshot` — the latter reads an
/// already-swapped back buffer and captures the previous frame.
///
/// Caller is responsible for having called `rl.initWindow` (and for closing
/// it afterward); this only owns the render texture and the per-shot export.
pub fn run(opts: Options, width: i32, height: i32, shots: []const Shot) !void {
    if (!opts.enabled) return;

    const rt = try rl.loadRenderTexture(width, height);
    defer rl.unloadRenderTexture(rt);

    var path_buf: [std.fs.max_path_bytes]u8 = undefined;
    for (shots) |shot| {
        rl.beginTextureMode(rt);
        shot.draw();
        rl.endTextureMode();

        var img = rl.loadImageFromTexture(rt.texture) catch continue;
        defer rl.unloadImage(img);
        rl.imageFlipVertical(&img); // render textures are y-flipped in GL
        const path = std.fmt.bufPrintZ(&path_buf, "{s}/{s}.png", .{ opts.out_dir, shot.name }) catch continue;
        _ = rl.exportImage(img, path);
    }
}

// ── Tests ────────────────────────────────────────────────────────────────
// `run` itself needs a live raylib window/GL context (see RULES.md — rl/ is
// compile-checked, not unit-tested), but the env-var parsing is pure std and
// gets covered here.

const testing = std.testing;

test "optionsFromEnviron: disabled when the enable var is unset or empty" {
    var map = std.process.Environ.Map.init(testing.allocator);
    defer map.deinit();

    var opts = optionsFromEnviron(&map, "APP_SHOT", "APP_SHOT_DIR");
    try testing.expect(!opts.enabled);
    try testing.expectEqualStrings("runs/shots", opts.out_dir);

    try map.put("APP_SHOT", "");
    opts = optionsFromEnviron(&map, "APP_SHOT", "APP_SHOT_DIR");
    try testing.expect(!opts.enabled);
}

test "optionsFromEnviron: enabled and custom dir from the environ map" {
    var map = std.process.Environ.Map.init(testing.allocator);
    defer map.deinit();
    try map.put("APP_SHOT", "1");
    try map.put("APP_SHOT_DIR", "/tmp/shots");

    const opts = optionsFromEnviron(&map, "APP_SHOT", "APP_SHOT_DIR");
    try testing.expect(opts.enabled);
    try testing.expectEqualStrings("/tmp/shots", opts.out_dir);
}
