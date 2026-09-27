const std = @import("std");
const contract = @import("../contract.zig");

/// Token spans are half-open byte offsets into the source slice the
/// `Lexer` was `init`ed with — zero-copy, no allocation. `tag(src)` slices
/// the text back out.
pub fn Token(comptime Tag: type) type {
    return struct {
        tag: Tag,
        start: usize,
        end: usize,

        pub fn text(tok: @This(), src: []const u8) []const u8 {
            return src[tok.start..tok.end];
        }
    };
}

/// Generic byte-driven lexer. `spec` is a tuple of `.{ "lexeme", Tag.value }`
/// pairs covering both keywords (matched against identifier-shaped tokens by
/// exact string equality) and punctuation/operators (matched greedily,
/// longest lexeme first, against raw source bytes). `Tag` is the caller's
/// token-tag enum and must declare `eof`, `ident`, `number`, `string` and
/// `invalid` members — the lexer reports those for anything not covered by
/// `spec`.
pub fn Lexer(comptime Tag: type, comptime spec: anytype) type {
    const Entry = struct { lexeme: []const u8, tag: Tag };

    const n = spec.len;
    const all: [n]Entry = blk: {
        var t: [n]Entry = undefined;
        for (spec, 0..) |e, i| t[i] = .{ .lexeme = e[0], .tag = e[1] };
        break :blk t;
    };

    const kw_count = blk: {
        var c: usize = 0;
        for (all) |e| {
            if (isIdentLexeme(e.lexeme)) c += 1;
        }
        break :blk c;
    };
    const keywords: [kw_count]Entry = blk: {
        var t: [kw_count]Entry = undefined;
        var i: usize = 0;
        for (all) |e| {
            if (isIdentLexeme(e.lexeme)) {
                t[i] = e;
                i += 1;
            }
        }
        break :blk t;
    };

    const punct_count = n - kw_count;
    const punct: [punct_count]Entry = blk: {
        var t: [punct_count]Entry = undefined;
        var i: usize = 0;
        for (all) |e| {
            if (!isIdentLexeme(e.lexeme)) {
                t[i] = e;
                i += 1;
            }
        }
        // Longest-first insertion sort: small comptime table, simplicity
        // over algorithmic cleverness.
        var j: usize = 1;
        while (j < punct_count) : (j += 1) {
            const key = t[j];
            var k = j;
            while (k > 0 and t[k - 1].lexeme.len < key.lexeme.len) : (k -= 1) {
                t[k] = t[k - 1];
            }
            t[k] = key;
        }
        break :blk t;
    };

    return struct {
        const Self = @This();
        pub const T = Token(Tag);

        src: []const u8,
        pos: usize = 0,

        pub fn init(src: []const u8) Self {
            return .{ .src = src };
        }

        /// Returns a token tagged `.eof` (zero-length, at `src.len`) once the
        /// input is exhausted; callers loop until they see it.
        pub fn next(lx: *Self) T {
            lx.skipTrivia();
            const start = lx.pos;
            if (lx.pos >= lx.src.len) return .{ .tag = .eof, .start = start, .end = start };

            const c = lx.src[lx.pos];

            if (isIdentStart(c)) {
                lx.pos += 1;
                while (lx.pos < lx.src.len and isIdentChar(lx.src[lx.pos])) lx.pos += 1;
                const word = lx.src[start..lx.pos];
                for (keywords) |kw| {
                    if (std.mem.eql(u8, kw.lexeme, word)) return .{ .tag = kw.tag, .start = start, .end = lx.pos };
                }
                return .{ .tag = .ident, .start = start, .end = lx.pos };
            }

            if (std.ascii.isDigit(c)) {
                lx.lexNumber();
                return .{ .tag = .number, .start = start, .end = lx.pos };
            }

            if (c == '"') {
                lx.lexString();
                return .{ .tag = .string, .start = start, .end = lx.pos };
            }

            inline for (punct) |p| {
                if (std.mem.startsWith(u8, lx.src[lx.pos..], p.lexeme)) {
                    lx.pos += p.lexeme.len;
                    return .{ .tag = p.tag, .start = start, .end = lx.pos };
                }
            }

            // Unrecognized byte: one-byte `.invalid` token; the cursor always advances.
            lx.pos += 1;
            return .{ .tag = .invalid, .start = start, .end = lx.pos };
        }

        /// Peek the next token without consuming it.
        pub fn peek(lx: *Self) T {
            var copy = lx.*;
            return copy.next();
        }

        fn skipTrivia(lx: *Self) void {
            while (lx.pos < lx.src.len) {
                const c = lx.src[lx.pos];
                if (c == ' ' or c == '\t' or c == '\n' or c == '\r') {
                    lx.pos += 1;
                } else if (c == '/' and lx.pos + 1 < lx.src.len and lx.src[lx.pos + 1] == '/') {
                    while (lx.pos < lx.src.len and lx.src[lx.pos] != '\n') lx.pos += 1;
                } else break;
            }
        }

        /// Decimal integers and floats only — no hex/octal/binary, no digit
        /// separators; callers wanting more can re-lex the `.number` text.
        fn lexNumber(lx: *Self) void {
            while (lx.pos < lx.src.len and std.ascii.isDigit(lx.src[lx.pos])) lx.pos += 1;
            if (lx.pos + 1 < lx.src.len and lx.src[lx.pos] == '.' and std.ascii.isDigit(lx.src[lx.pos + 1])) {
                lx.pos += 1;
                while (lx.pos < lx.src.len and std.ascii.isDigit(lx.src[lx.pos])) lx.pos += 1;
            }
            if (lx.pos < lx.src.len and (lx.src[lx.pos] == 'e' or lx.src[lx.pos] == 'E')) {
                var k = lx.pos + 1;
                if (k < lx.src.len and (lx.src[k] == '+' or lx.src[k] == '-')) k += 1;
                if (k < lx.src.len and std.ascii.isDigit(lx.src[k])) {
                    lx.pos = k;
                    while (lx.pos < lx.src.len and std.ascii.isDigit(lx.src[lx.pos])) lx.pos += 1;
                }
            }
        }

        /// An unterminated string runs to end-of-input rather than looping
        /// forever; the span still includes the opening quote so callers can
        /// detect the missing closer by checking the last byte.
        fn lexString(lx: *Self) void {
            lx.pos += 1; // opening quote
            while (lx.pos < lx.src.len) {
                const c = lx.src[lx.pos];
                if (c == '"') {
                    lx.pos += 1;
                    return;
                }
                if (c == '\\' and lx.pos + 1 < lx.src.len) {
                    lx.pos += 2;
                } else {
                    lx.pos += 1;
                }
            }
        }

        fn isIdentStart(c: u8) bool {
            return std.ascii.isAlphabetic(c) or c == '_';
        }
        fn isIdentChar(c: u8) bool {
            return std.ascii.isAlphanumeric(c) or c == '_';
        }
    };
}

fn isIdentLexeme(lexeme: []const u8) bool {
    contract.require(lexeme.len > 0, "lex: spec entry has empty lexeme");
    for (lexeme, 0..) |c, i| {
        const ok = std.ascii.isAlphanumeric(c) or c == '_';
        if (!ok) return false;
        if (i == 0 and std.ascii.isDigit(c)) return false;
    }
    return true;
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

const TestTag = enum {
    eof,
    ident,
    number,
    string,
    invalid,
    kw_if,
    kw_else,
    kw_let,
    plus,
    eq,
    eq_eq,
    eq_eq_eq,
    lt,
    lt_eq,
    lparen,
    rparen,
    lbrace,
    rbrace,
    semi,
};

const test_spec = .{
    .{ "if", TestTag.kw_if },
    .{ "else", TestTag.kw_else },
    .{ "let", TestTag.kw_let },
    .{ "+", TestTag.plus },
    .{ "=", TestTag.eq },
    .{ "==", TestTag.eq_eq },
    .{ "===", TestTag.eq_eq_eq },
    .{ "<", TestTag.lt },
    .{ "<=", TestTag.lt_eq },
    .{ "(", TestTag.lparen },
    .{ ")", TestTag.rparen },
    .{ "{", TestTag.lbrace },
    .{ "}", TestTag.rbrace },
    .{ ";", TestTag.semi },
};

const TestLexer = Lexer(TestTag, test_spec);

fn collect(src: []const u8, out: []TestLexer.T) []TestLexer.T {
    var lx = TestLexer.init(src);
    var i: usize = 0;
    while (true) {
        const t = lx.next();
        out[i] = t;
        i += 1;
        if (t.tag == .eof) break;
    }
    return out[0..i];
}

test Lexer {
    const src =
        \\let x = 1.5e2; // trailing comment
        \\if x <= 3 {
        \\  x === "a\"b"
        \\} else {
        \\  x + 2
        \\}
    ;
    var buf: [64]TestLexer.T = undefined;
    const toks = collect(src, &buf);

    const want = [_]TestTag{
        .kw_let, .ident,    .eq,     .number, .semi,
        .kw_if,  .ident,    .lt_eq,  .number, .lbrace,
        .ident,  .eq_eq_eq, .string, .rbrace, .kw_else,
        .lbrace, .ident,    .plus,   .number, .rbrace,
        .eof,
    };
    try testing.expectEqual(want.len, toks.len);
    for (toks, want) |t, w| try testing.expectEqual(w, t.tag);

    // Spot-check exact spans for a few tokens.
    const let_tok = toks[0];
    try testing.expectEqualStrings("let", let_tok.text(src));

    const num_tok = toks[3];
    try testing.expectEqualStrings("1.5e2", num_tok.text(src));

    const str_tok = toks[12];
    try testing.expectEqualStrings("\"a\\\"b\"", str_tok.text(src));
}

test "longest-match punctuation: === vs == vs =, <= vs <" {
    var lx = TestLexer.init("=== == = <= <");
    try testing.expectEqualStrings("===", lx.next().text("=== == = <= <"));
    try testing.expectEqualStrings("==", lx.next().text("=== == = <= <"));
    try testing.expectEqualStrings("=", lx.next().text("=== == = <= <"));
    try testing.expectEqualStrings("<=", lx.next().text("=== == = <= <"));
    try testing.expectEqualStrings("<", lx.next().text("=== == = <= <"));
    try testing.expectEqual(TestTag.eof, lx.next().tag);
}

test "peek does not consume" {
    var lx = TestLexer.init("if x");
    const p1 = lx.peek();
    const p2 = lx.peek();
    try testing.expectEqual(TestTag.kw_if, p1.tag);
    try testing.expectEqual(p1.tag, p2.tag);
    try testing.expectEqual(p1.start, p2.start);
    const n1 = lx.next();
    try testing.expectEqual(TestTag.kw_if, n1.tag);
    try testing.expectEqual(TestTag.ident, lx.next().tag);
}

test "Token.text slices the source by span" {
    const src = "let answer = 42;";
    var lx = TestLexer.init(src);
    const kw = lx.next();
    const ident = lx.next();
    try testing.expectEqualStrings("let", kw.text(src));
    try testing.expectEqualStrings("answer", ident.text(src));
}

test "keyword vs identifier: exact match after ident scan" {
    var lx = TestLexer.init("iffy if_ ifelse if");
    try testing.expectEqual(TestTag.ident, lx.next().tag);
    try testing.expectEqual(TestTag.ident, lx.next().tag);
    try testing.expectEqual(TestTag.ident, lx.next().tag);
    try testing.expectEqual(TestTag.kw_if, lx.next().tag);
}

test "numbers: int, float, exponent" {
    var lx = TestLexer.init("42 3.14 2e10 5e-3");
    for ([_][]const u8{ "42", "3.14", "2e10", "5e-3" }) |want| {
        const t = lx.next();
        try testing.expectEqual(TestTag.number, t.tag);
        try testing.expectEqualStrings(want, t.text("42 3.14 2e10 5e-3"));
    }
}

test "string with escapes, unterminated string reaches eof safely" {
    var lx = TestLexer.init("\"hello\\nworld\" \"unterminated");
    const t1 = lx.next();
    try testing.expectEqual(TestTag.string, t1.tag);
    try testing.expectEqualStrings("\"hello\\nworld\"", t1.text("\"hello\\nworld\" \"unterminated"));
    const t2 = lx.next();
    try testing.expectEqual(TestTag.string, t2.tag);
    try testing.expectEqual(@as(usize, 28), t2.end);
    try testing.expectEqual(TestTag.eof, lx.next().tag);
}

test "line comments are skipped as trivia" {
    var lx = TestLexer.init("x // comment ignored\ny");
    try testing.expectEqual(TestTag.ident, lx.next().tag);
    const y = lx.next();
    try testing.expectEqual(TestTag.ident, y.tag);
    try testing.expectEqualStrings("y", y.text("x // comment ignored\ny"));
    try testing.expectEqual(TestTag.eof, lx.next().tag);
}

test "whitespace-only and empty source yield only eof" {
    var lx = TestLexer.init("   \n\t  ");
    try testing.expectEqual(TestTag.eof, lx.next().tag);
    var lx2 = TestLexer.init("");
    try testing.expectEqual(TestTag.eof, lx2.next().tag);
}

test "unrecognized byte yields one-byte invalid token" {
    const L = Lexer(TestTag, .{.{ "+", .plus }});
    var lx = L.init("a @ b");
    try testing.expectEqual(TestTag.ident, lx.next().tag);
    const bad = lx.next();
    try testing.expectEqual(TestTag.invalid, bad.tag);
    try testing.expectEqual(@as(usize, 1), bad.end - bad.start);
    try testing.expectEqual(TestTag.ident, lx.next().tag);
    try testing.expectEqual(TestTag.eof, lx.next().tag);
}
