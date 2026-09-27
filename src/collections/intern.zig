//! Interner — a string interner that hands out small `Handle`s in place of
//! `[]const u8` so equal strings compare equal in O(1) and callers can stash
//! a handle instead of carrying a slice (and its lifetime) around. Every
//! interned string's bytes are concatenated into one arena and the map is
//! keyed by offset into it, so one allocation is amortized over the whole
//! arena rather than paying one allocation per interned string.
//!
//! Arena growth can `realloc` and move the backing buffer, so a `[]const u8`
//! returned by `get` is only valid until the next `intern` call — do not
//! cache it across further interning. `Handle` values (arena offsets, not
//! pointers) stay valid for the Interner's whole lifetime.

const std = @import("std");
const contract = @import("../contract.zig");

/// Opaque handle to an interned string. Backed by its byte offset into the
/// arena, but non-exhaustive and integer-tagged so callers can't do offset
/// arithmetic on it directly or assume adjacency between handles.
pub const Handle = enum(u32) { _ };

const Entry = struct { off: u32, len: u32 };

/// Hashes/compares a `u32` key already in the map by resolving it to bytes
/// in the arena; used to rehash existing keys on grow.
const IndexContext = struct {
    entries: *const std.ArrayList(Entry),
    arena: *const std.ArrayList(u8),

    fn bytesOf(self: @This(), key: u32) []const u8 {
        const e = self.entries.items[key];
        return self.arena.items[e.off..][0..e.len];
    }

    pub fn hash(self: @This(), key: u32) u64 {
        return std.hash.Wyhash.hash(0, self.bytesOf(key));
    }

    pub fn eql(self: @This(), a: u32, b: u32) bool {
        return std.mem.eql(u8, self.bytesOf(a), self.bytesOf(b));
    }
};

/// Hashes/compares a probe `[]const u8` against a stored `u32` key resolved
/// through the arena, so lookup never has to allocate/copy the probe string.
const SliceAdapter = struct {
    entries: *const std.ArrayList(Entry),
    arena: *const std.ArrayList(u8),

    fn bytesOf(self: @This(), key: u32) []const u8 {
        const e = self.entries.items[key];
        return self.arena.items[e.off..][0..e.len];
    }

    pub fn hash(_: @This(), probe: []const u8) u64 {
        return std.hash.Wyhash.hash(0, probe);
    }

    pub fn eql(self: @This(), probe: []const u8, key: u32) bool {
        return std.mem.eql(u8, probe, self.bytesOf(key));
    }
};

/// String interner: dedups byte strings into handles comparable with `==`.
pub const Interner = struct {
    /// Concatenated bytes of every interned string.
    arena: std.ArrayList(u8) = .empty,
    /// `Handle` -> span into `arena`; index is `@intFromEnum(handle)`.
    entries: std.ArrayList(Entry) = .empty,
    /// Handle index deduped by the bytes it denotes (see `IndexContext`).
    map: std.HashMapUnmanaged(u32, void, IndexContext, std.hash_map.default_max_load_percentage) = .empty,

    pub const empty: Interner = .{};

    pub fn deinit(self: *Interner, gpa: std.mem.Allocator) void {
        self.arena.deinit(gpa);
        self.entries.deinit(gpa);
        self.map.deinit(gpa);
        self.* = undefined;
    }

    fn indexContext(self: *const Interner) IndexContext {
        return .{ .entries = &self.entries, .arena = &self.arena };
    }

    fn sliceAdapter(self: *const Interner) SliceAdapter {
        return .{ .entries = &self.entries, .arena = &self.arena };
    }

    /// Look up `s` without interning it. Returns `null` if `s` has not been
    /// interned yet.
    pub fn lookup(self: *const Interner, s: []const u8) ?Handle {
        const idx = self.map.getKeyAdapted(s, self.sliceAdapter()) orelse return null;
        return @enumFromInt(idx);
    }

    /// Intern `s`, returning its existing `Handle` if `s` was already
    /// interned, or appending it to the arena and returning a fresh `Handle`
    /// otherwise. `s` is never stored more than once.
    pub fn intern(self: *Interner, gpa: std.mem.Allocator, s: []const u8) !Handle {
        contract.require(self.entries.items.len < std.math.maxInt(u32), "intern: interner full");

        const gop = try self.map.getOrPutContextAdapted(
            gpa,
            s,
            self.sliceAdapter(),
            self.indexContext(),
        );
        if (gop.found_existing) return @enumFromInt(gop.key_ptr.*);
        errdefer _ = self.map.removeContext(gop.key_ptr.*, self.indexContext());

        const off: u32 = @intCast(self.arena.items.len);
        try self.arena.appendSlice(gpa, s);
        errdefer self.arena.shrinkRetainingCapacity(off);

        const idx: u32 = @intCast(self.entries.items.len);
        try self.entries.append(gpa, .{ .off = off, .len = @intCast(s.len) });
        gop.key_ptr.* = idx;
        return @enumFromInt(idx);
    }

    /// Resolve `h` back to the bytes it denotes. Zero-copy: the slice points
    /// into the arena and is valid only until the next `intern` call.
    pub fn get(self: *const Interner, h: Handle) []const u8 {
        const idx = @intFromEnum(h);
        contract.require(idx < self.entries.items.len, "get: Handle out of range");
        const e = self.entries.items[idx];
        return self.arena.items[e.off..][0..e.len];
    }
};

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test Interner {
    var itn: Interner = .empty;
    defer itn.deinit(testing.allocator);

    const a = try itn.intern(testing.allocator, "hello");
    const len_after_first = itn.arena.items.len;
    const b = try itn.intern(testing.allocator, "hello");

    try testing.expectEqual(a, b);
    try testing.expectEqual(len_after_first, itn.arena.items.len);
    try testing.expectEqualStrings("hello", itn.get(a));
}

test "interning distinct strings yields distinct handles that resolve correctly" {
    var itn: Interner = .empty;
    defer itn.deinit(testing.allocator);

    const a = try itn.intern(testing.allocator, "foo");
    const b = try itn.intern(testing.allocator, "bar");
    const c = try itn.intern(testing.allocator, "foobar");

    try testing.expect(a != b);
    try testing.expect(b != c);
    try testing.expect(a != c);
    try testing.expectEqualStrings("foo", itn.get(a));
    try testing.expectEqualStrings("bar", itn.get(b));
    try testing.expectEqualStrings("foobar", itn.get(c));
}

test "lookup finds a string only after it has been interned" {
    var itn: Interner = .empty;
    defer itn.deinit(testing.allocator);

    try testing.expectEqual(@as(?Handle, null), itn.lookup("not yet interned"));

    const h = try itn.intern(testing.allocator, "not yet interned");
    try testing.expectEqual(h, itn.lookup("not yet interned").?);
    try testing.expectEqual(@as(?Handle, null), itn.lookup("still absent"));
}

test "deinit frees a handful of interned strings with no leaks" {
    var itn: Interner = .empty;
    defer itn.deinit(testing.allocator);

    _ = try itn.intern(testing.allocator, "alpha");
    _ = try itn.intern(testing.allocator, "beta");
    _ = try itn.intern(testing.allocator, "gamma");
    _ = try itn.intern(testing.allocator, "alpha"); // dedup path
    _ = try itn.intern(testing.allocator, "delta");
}

test "empty string interns and resolves" {
    var itn: Interner = .empty;
    defer itn.deinit(testing.allocator);

    const h = try itn.intern(testing.allocator, "");
    try testing.expectEqualStrings("", itn.get(h));
    try testing.expectEqual(h, (try itn.intern(testing.allocator, "")));
}
