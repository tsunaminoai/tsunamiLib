const std = @import("std");
const contract = @import("../contract.zig");

/// dlopen-based hot reload of a single symbol, polled by mtime.
///
/// `Fn` is the function-pointer type looked up in the shared library (e.g.
/// `fn (*State) callconv(.c) void`); `sym` is its comptime-known,
/// nul-terminated symbol name. Callers own the `std.Io` instance and pass it
/// explicitly to every call that needs to touch the filesystem.
pub fn HotReload(comptime Fn: type, comptime sym: [:0]const u8) type {
    return struct {
        const Self = @This();

        lib: ?std.DynLib = null,
        fn_ptr: ?Fn = null,
        last_mtime: i128 = 0,
        gen: u32 = 0,
        shadow_buf: [std.fs.max_path_bytes]u8 = undefined,
        shadow_len: usize = 0,

        pub const Error = std.DynLib.Error || error{ SymbolNotFound, CopyFailed, NameTooLong };

        /// Loads `path` for the first time, populating `fn_ptr` and
        /// `last_mtime`. Leaves the instance untouched (no lib open, no
        /// function pointer) on failure so the caller can fall back to a
        /// stub and retry later via `poll`.
        pub fn init(io: std.Io, path: []const u8) Error!Self {
            contract.require(path.len > 0, "hotreload.init: empty path");
            var self: Self = .{};
            try self.load(io, path);
            self.last_mtime = statMtime(io, path) orelse 0;
            return self;
        }

        pub fn deinit(self: *Self, io: std.Io) void {
            if (self.lib) |*lib| lib.close();
            self.dropShadow(io);
            self.* = undefined;
        }

        /// Current function pointer, or `null` if nothing has loaded
        /// successfully yet.
        pub fn get(self: *const Self) ?Fn {
            return self.fn_ptr;
        }

        /// Stats `path`; if its mtime advanced since the last successful
        /// load (or `init`), reloads it. Returns whether a reload happened.
        /// A stat failure (file missing mid-rebuild, etc.) or a dlopen/
        /// lookup failure is swallowed: this is an expected, transient state
        /// while a background build is in progress, so `poll` just leaves
        /// the previous `fn_ptr` in place and reports no reload.
        pub fn poll(self: *Self, io: std.Io, path: []const u8) bool {
            contract.require(path.len > 0, "hotreload.poll: empty path");
            const mtime = statMtime(io, path) orelse return false;
            if (mtime <= self.last_mtime) return false;
            self.load(io, path) catch return false;
            self.last_mtime = mtime;
            return true;
        }

        /// dlopen caches handles by path, so reopening the same path while the
        /// old handle is live returns the OLD code. Each load copies the
        /// library to a fresh `<path>.hot<gen>` shadow and opens that instead.
        fn load(self: *Self, io: std.Io, path: []const u8) Error!void {
            var name_buf: [std.fs.max_path_bytes]u8 = undefined;
            const shadow = std.fmt.bufPrint(&name_buf, "{s}.hot{d}", .{ path, self.gen }) catch return error.NameTooLong;
            const cwd = std.Io.Dir.cwd();
            cwd.copyFile(path, cwd, shadow, io, .{}) catch |e| return if (e == error.FileNotFound) error.FileNotFound else error.CopyFailed;
            var new_lib = std.DynLib.open(shadow) catch |e| {
                cwd.deleteFile(io, shadow) catch {};
                return e;
            };
            const looked_up = new_lib.lookup(Fn, sym) orelse {
                new_lib.close();
                cwd.deleteFile(io, shadow) catch {};
                return error.SymbolNotFound;
            };
            if (self.lib) |*old| old.close();
            self.dropShadow(io);
            @memcpy(self.shadow_buf[0..shadow.len], shadow);
            self.shadow_len = shadow.len;
            self.gen +%= 1;
            self.lib = new_lib;
            self.fn_ptr = looked_up;
        }

        fn dropShadow(self: *Self, io: std.Io) void {
            if (self.shadow_len == 0) return;
            std.Io.Dir.cwd().deleteFile(io, self.shadow_buf[0..self.shadow_len]) catch {};
            self.shadow_len = 0;
        }

        fn statMtime(io: std.Io, path: []const u8) ?i128 {
            const stat = std.Io.Dir.cwd().statFile(io, path, .{}) catch return null;
            return stat.mtime.nanoseconds;
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

/// `statMtime`, `poll`, and `init`'s mtime path all key off `Dir.cwd()` +
/// `sub_path`, so tests drive them through an absolute path inside a tmp dir
/// rather than relying on the test process's actual working directory.
fn absPath(dir: std.Io.Dir, io: std.Io, buf: []u8, sub_path: []const u8) ![]const u8 {
    var dir_buf: [std.fs.max_path_bytes]u8 = undefined;
    const dir_len = try dir.realPath(io, &dir_buf);
    return std.fmt.bufPrint(buf, "{s}/{s}", .{ dir_buf[0..dir_len], sub_path });
}

test "statMtime sees a freshly created file" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(testing.io, .{ .sub_path = "lib.so", .data = "v1" });

    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try absPath(tmp.dir, testing.io, &buf, "lib.so");

    const H = HotReload(*const fn () callconv(.c) void, "run");
    const mtime = H.statMtime(testing.io, path) orelse return error.TestUnexpectedResult;
    try testing.expect(mtime > 0);
}

test "statMtime returns null for a missing file" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();

    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try absPath(tmp.dir, testing.io, &buf, "does-not-exist.so");

    const H = HotReload(*const fn () callconv(.c) void, "run");
    try testing.expect(H.statMtime(testing.io, path) == null);
}

test "poll detects a bumped mtime without a real shared library" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();

    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try absPath(tmp.dir, testing.io, &buf, "watched.so");

    try tmp.dir.writeFile(testing.io, .{ .sub_path = "watched.so", .data = "v1" });

    const H = HotReload(*const fn () callconv(.c) void, "run");
    const first = H.statMtime(testing.io, path) orelse return error.TestUnexpectedResult;

    // Force a distinct, later mtime: on filesystems with coarse mtime
    // resolution a same-tick rewrite can appear unchanged.
    std.Io.sleep(testing.io, .fromMilliseconds(10), .real) catch {};
    try tmp.dir.writeFile(testing.io, .{ .sub_path = "watched.so", .data = "v2, longer payload" });

    const second = H.statMtime(testing.io, path) orelse return error.TestUnexpectedResult;
    try testing.expect(second > first);

    // Now drive the actual `poll` entry point (still without a real
    // library, so the load attempt fails and poll reports no reload — but
    // last_mtime must not be silently bumped either).
    var hr: H = .{ .last_mtime = first };
    const reloaded = hr.poll(testing.io, path);
    try testing.expect(!reloaded);
    try testing.expectEqual(first, hr.last_mtime);
}

test "poll is a no-op when mtime is unchanged" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();

    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try absPath(tmp.dir, testing.io, &buf, "stable.so");
    try tmp.dir.writeFile(testing.io, .{ .sub_path = "stable.so", .data = "v1" });

    const H = HotReload(*const fn () callconv(.c) void, "run");
    const mtime = H.statMtime(testing.io, path) orelse return error.TestUnexpectedResult;

    var hr: H = .{ .last_mtime = mtime };
    try testing.expect(!hr.poll(testing.io, path));
    try testing.expectEqual(mtime, hr.last_mtime);
}

test "poll and init report no reload / propagate error for a missing file" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();

    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try absPath(tmp.dir, testing.io, &buf, "missing.so");

    const H = HotReload(*const fn () callconv(.c) void, "run");
    var hr: H = .{};
    try testing.expect(!hr.poll(testing.io, path));
    try testing.expectEqual(@as(i128, 0), hr.last_mtime);

    try testing.expectError(error.FileNotFound, H.init(testing.io, path));
}

test "get returns null before any successful load" {
    const H = HotReload(*const fn () callconv(.c) void, "run");
    var hr: H = .{};
    defer hr.deinit(testing.io);
    try testing.expect(hr.get() == null);
}

test "failed load removes its shadow copy" {
    var tmp = testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(testing.io, .{ .sub_path = "notalib.so", .data = "not an ELF" });
    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const path = try absPath(tmp.dir, testing.io, &buf, "notalib.so");
    const H = HotReload(*const fn () callconv(.c) void, "run");
    try testing.expect(std.meta.isError(H.init(testing.io, path)));
    try testing.expectError(error.FileNotFound, tmp.dir.statFile(testing.io, "notalib.so.hot0", .{}));
}
