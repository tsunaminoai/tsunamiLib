const std = @import("std");
const rl = @import("raylib");

/// Every drawable color in a UI comes from one of these fields, so an entire
/// app re-tints by swapping which `Palette` is active.
pub const Palette = struct {
    background: rl.Color,
    panel: rl.Color,
    text: rl.Color,
    text_soft: rl.Color,
    accent: rl.Color,
    button: rl.Color,
    button_hover: rl.Color,
    disabled: rl.Color,
};

pub const dark: Palette = .{
    .background = .{ .r = 18, .g = 20, .b = 24, .a = 255 },
    .panel = .{ .r = 30, .g = 33, .b = 38, .a = 235 },
    .text = .{ .r = 236, .g = 236, .b = 232, .a = 255 },
    .text_soft = .{ .r = 168, .g = 168, .b = 164, .a = 255 },
    .accent = .{ .r = 224, .g = 178, .b = 74, .a = 255 },
    .button = .{ .r = 52, .g = 58, .b = 68, .a = 255 },
    .button_hover = .{ .r = 74, .g = 82, .b = 94, .a = 255 },
    .disabled = .{ .r = 70, .g = 72, .b = 76, .a = 255 },
};

pub const light: Palette = .{
    .background = .{ .r = 244, .g = 244, .b = 240, .a = 255 },
    .panel = .{ .r = 255, .g = 255, .b = 255, .a = 235 },
    .text = .{ .r = 24, .g = 24, .b = 24, .a = 255 },
    .text_soft = .{ .r = 96, .g = 96, .b = 96, .a = 255 },
    .accent = .{ .r = 32, .g = 110, .b = 200, .a = 255 },
    .button = .{ .r = 224, .g = 224, .b = 224, .a = 255 },
    .button_hover = .{ .r = 200, .g = 200, .b = 200, .a = 255 },
    .disabled = .{ .r = 210, .g = 210, .b = 210, .a = 255 },
};

/// A complete, swappable theme: palette plus render-affecting handles.
pub const Theme = struct {
    palette: Palette,
    font: ?rl.Font = null,
    font_spacing: f32 = 1.0,
};

pub const dark_theme: Theme = .{ .palette = dark };
pub const light_theme: Theme = .{ .palette = light };

/// The active theme. Callers read `current.palette.*` to draw; `apply`
/// switches it (e.g. at the top of a frame or on a settings toggle).
pub var current: Theme = dark_theme;

pub fn apply(t: Theme) void {
    current = t;
}

/// A color with its alpha replaced (0..1, clamped), independent of `apply`.
pub fn withAlpha(c: rl.Color, alpha: f32) rl.Color {
    const a = std.math.clamp(alpha, 0.0, 1.0);
    return .{ .r = c.r, .g = c.g, .b = c.b, .a = @intFromFloat(a * 255.0) };
}

/// A rounded, labeled button using the active palette, styled by `hovered`
/// and `enabled`.
pub fn drawButton(rect: rl.Rectangle, label: [:0]const u8, hovered: bool, enabled: bool) void {
    const p = current.palette;
    const bg = if (!enabled) p.disabled else if (hovered) p.button_hover else p.button;
    rl.drawRectangleRounded(rect, 0.3, 6, bg);
    rl.drawRectangleRoundedLinesEx(rect, 0.3, 6, 2.0, if (enabled) p.accent else p.text_soft);
    const size: i32 = @intFromFloat(rect.height * 0.42);
    const tw = rl.measureText(label, size);
    rl.drawText(
        label,
        @as(i32, @intFromFloat(rect.x + rect.width / 2.0)) - @divTrunc(tw, 2),
        @intFromFloat(rect.y + (rect.height - @as(f32, @floatFromInt(size))) / 2.0),
        size,
        if (enabled) p.text else p.text_soft,
    );
}

// No tests: drawing calls need a live raylib window/GL context. See
// RULES.md — rl/ is compile-checked, not unit-tested.
