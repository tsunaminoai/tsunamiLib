const std = @import("std");

/// UTF-8 aware typewriter text reveal: advances a fractional character count
/// each `update(dt)` and exposes the byte offset to slice the source text at
/// for drawing. Pure logic — no rendering, no allocation.
pub const Typewriter = struct {
    /// Identity of the text being revealed (ptr+len), so switching to a new
    /// string resets the reveal instead of carrying over stale progress.
    text_ptr: usize = 0,
    text_len: usize = 0,
    /// Codepoints revealed so far, fractional so partial-second updates
    /// accumulate correctly.
    chars: f32 = 0,
    chars_per_second: f32 = 30,

    /// Advance the reveal toward the end of `text`, resetting to the start
    /// whenever `text` identifies a different string than last call.
    pub fn update(self: *Typewriter, text: []const u8, dt: f32) void {
        const id = @intFromPtr(text.ptr);
        if (id != self.text_ptr or text.len != self.text_len) {
            self.text_ptr = id;
            self.text_len = text.len;
            self.chars = 0;
        }
        const total: f32 = @floatFromInt(countCodepoints(text));
        self.chars = @min(self.chars + dt * self.chars_per_second, total);
    }

    /// The prefix of `text` revealed so far, cut on a codepoint boundary.
    /// `text` must be the same string most recently passed to `update`.
    pub fn revealed(self: *const Typewriter, text: []const u8) []const u8 {
        const n: usize = @intFromFloat(self.chars);
        return text[0..byteOffsetForCodepoints(text, n)];
    }

    /// Whether every codepoint of `text` has been revealed.
    pub fn complete(self: *const Typewriter, text: []const u8) bool {
        return self.revealed(text).len >= text.len;
    }

    /// Jump straight to fully revealed (e.g. on a skip click).
    pub fn finish(self: *Typewriter, text: []const u8) void {
        self.chars = @floatFromInt(countCodepoints(text));
    }
};

fn countCodepoints(text: []const u8) usize {
    var n: usize = 0;
    var i: usize = 0;
    while (i < text.len) {
        i += std.unicode.utf8ByteSequenceLength(text[i]) catch 1;
        n += 1;
    }
    return n;
}

fn byteOffsetForCodepoints(text: []const u8, n: usize) usize {
    var count: usize = 0;
    var i: usize = 0;
    while (i < text.len and count < n) {
        i += std.unicode.utf8ByteSequenceLength(text[i]) catch 1;
        count += 1;
    }
    return i;
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "reveals incrementally and completes" {
    var tw: Typewriter = .{ .chars_per_second = 10 };
    const text = "hello world";
    tw.update(text, 0.5); // 5 chars
    try testing.expectEqualStrings("hello", tw.revealed(text));
    try testing.expect(!tw.complete(text));
    tw.update(text, 10.0); // way past the end, clamps
    try testing.expectEqualStrings(text, tw.revealed(text));
    try testing.expect(tw.complete(text));
}

test "finish jumps straight to the end" {
    var tw: Typewriter = .{};
    const text = "abcdef";
    tw.finish(text);
    try testing.expectEqualStrings(text, tw.revealed(text));
    try testing.expect(tw.complete(text));
}

test "switching text identity resets progress" {
    var tw: Typewriter = .{ .chars_per_second = 100 };
    const a = "first string";
    tw.update(a, 1.0);
    try testing.expect(tw.complete(a));
    const b = "a different second string";
    tw.update(b, 0.0);
    try testing.expectEqualStrings("", tw.revealed(b));
}

test "UTF-8 aware: reveal cuts on codepoint boundaries" {
    var tw: Typewriter = .{ .chars_per_second = 10 };
    const text = "h\u{00e9}llo"; // 'é' is 2 bytes, 5 codepoints total
    tw.update(text, 0.3); // 3 codepoints: h, é, l
    const rev = tw.revealed(text);
    try testing.expectEqual(@as(usize, 4), rev.len); // 1 + 2 + 1 bytes
    try testing.expect(std.unicode.utf8ValidateSlice(rev));
    tw.update(text, 10.0);
    try testing.expectEqualStrings(text, tw.revealed(text));
}
