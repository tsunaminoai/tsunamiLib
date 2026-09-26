//! Best-effort persistence of a byte blob: a small file under a caller-given
//! directory natively, `localStorage` under Emscripten. Every failure comes
//! back as an error (a missing file is `error.NotFound`) rather than a panic.

const std = @import("std");
const builtin = @import("builtin");
const contract = @import("../contract.zig");

const is_web = builtin.target.os.tag == .emscripten;

/// Errors `save`/`load` can return on the native backend. Deliberately small:
/// callers branch on `error.NotFound` and otherwise treat any other member as
/// "couldn't persist this time."
pub const Error = error{
    /// No save data under this name (first run, or it was never written).
    NotFound,
    /// The caller's buffer was too small to hold the saved bytes.
    BufferTooSmall,
    /// Any other I/O failure (permissions, disk full, etc).
    Failed,
};

// ── web: localStorage via emscripten_run_script ──────────────────────────

extern fn emscripten_run_script(script: [*:0]const u8) void;
extern fn emscripten_run_script_string(script: [*:0]const u8) [*:0]const u8;

/// Hex, because the script bridge only moves C strings and a save blob is
/// binary. 2 chars per byte.
fn toHex(bytes: []const u8, out: []u8) []const u8 {
    const digits = "0123456789abcdef";
    var n: usize = 0;
    for (bytes) |b| {
        if (n + 2 > out.len) break;
        out[n] = digits[b >> 4];
        out[n + 1] = digits[b & 0xF];
        n += 2;
    }
    return out[0..n];
}

fn fromHex(hex: []const u8, out: []u8) []const u8 {
    var n: usize = 0;
    var i: usize = 0;
    while (i + 1 < hex.len and n < out.len) : (i += 2) {
        const hi = std.fmt.charToDigit(hex[i], 16) catch break;
        const lo = std.fmt.charToDigit(hex[i + 1], 16) catch break;
        out[n] = (@as(u8, hi) << 4) | lo;
        n += 1;
    }
    return out[0..n];
}

// ── public API ─────────────────────────────────────────────────────────────

/// Persist `bytes` under `name` in `dir` (native) or under `name` as a
/// localStorage key (web, `dir` is ignored there). Best-effort: any I/O
/// failure is reported as `error.Failed` rather than surfaced in detail.
pub fn save(gpa: std.mem.Allocator, io: std.Io, dir: std.Io.Dir, name: []const u8, bytes: []const u8) Error!void {
    contract.require(name.len > 0, "save.save: empty name");

    if (comptime is_web) {
        const hex = gpa.alloc(u8, bytes.len * 2) catch return error.Failed;
        defer gpa.free(hex);
        const h = toHex(bytes, hex);

        const script = std.fmt.allocPrintSentinel(
            gpa,
            "try{{localStorage.setItem('{s}','{s}')}}catch(e){{}}",
            .{ name, h },
            0,
        ) catch return error.Failed;
        defer gpa.free(script);
        emscripten_run_script(script.ptr);
    } else {
        dir.writeFile(io, .{ .sub_path = name, .data = bytes }) catch return error.Failed;
    }
}

/// Load previously-saved bytes for `name` into `out`, returning the
/// used prefix. `error.NotFound` means there is nothing saved under this
/// name (a fresh run, not a failure); `error.BufferTooSmall` means `out`
/// couldn't hold the saved data.
pub fn load(gpa: std.mem.Allocator, io: std.Io, dir: std.Io.Dir, name: []const u8, out: []u8) Error![]u8 {
    contract.require(name.len > 0, "save.load: empty name");

    if (comptime is_web) {
        const key_script = std.fmt.allocPrintSentinel(
            gpa,
            "(function(){{try{{return localStorage.getItem('{s}')||''}}catch(e){{return ''}}}})()",
            .{name},
            0,
        ) catch return error.Failed;
        defer gpa.free(key_script);

        const raw = emscripten_run_script_string(key_script.ptr);
        const hex = std.mem.span(raw);
        if (hex.len == 0) return error.NotFound;

        const decoded = fromHex(hex, out);
        if (decoded.len == out.len and hex.len > out.len * 2) return error.BufferTooSmall;
        return decoded;
    } else {
        return dir.readFile(io, name, out) catch |err| switch (err) {
            error.FileNotFound => error.NotFound,
            else => error.Failed,
        };
    }
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "save then load round-trips the same bytes" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();

    const payload = "hello, save file";
    try save(testing.allocator, testing.io, tmp.dir, "slot", payload);

    var buf: [64]u8 = undefined;
    const back = try load(testing.allocator, testing.io, tmp.dir, "slot", &buf);
    try testing.expectEqualSlices(u8, payload, back);
}

test "load from a nonexistent file returns error.NotFound" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();

    var buf: [64]u8 = undefined;
    const result = load(testing.allocator, testing.io, tmp.dir, "does-not-exist", &buf);
    try testing.expectError(error.NotFound, result);
}

test "save overwrites a previous save under the same name" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();

    try save(testing.allocator, testing.io, tmp.dir, "slot", "first");
    try save(testing.allocator, testing.io, tmp.dir, "slot", "second value");

    var buf: [64]u8 = undefined;
    const back = try load(testing.allocator, testing.io, tmp.dir, "slot", &buf);
    try testing.expectEqualSlices(u8, "second value", back);
}
