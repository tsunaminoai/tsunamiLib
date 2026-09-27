//! A tiny arithmetic calculator: `text.lex.Lexer` tokenizes, `text.pratt.Parser`
//! drives precedence climbing, and a Builder evaluates directly to f64.
const std = @import("std");
const ts = @import("tsunami");

const Tag = enum { eof, ident, number, string, invalid, plus, minus, star, slash, caret, lparen, rparen };

const Lexer = ts.text.lex.Lexer(Tag, .{
    .{ "+", Tag.plus },
    .{ "-", Tag.minus },
    .{ "*", Tag.star },
    .{ "/", Tag.slash },
    .{ "^", Tag.caret },
    .{ "(", Tag.lparen },
    .{ ")", Tag.rparen },
});
const Tok = Lexer.T;

const Parser = ts.text.pratt.Parser(Tag, .{
    .prefix = &.{.{ .minus, 7 }},
    .infix = &.{
        .{ .plus, 1, 1 },
        .{ .minus, 1, 1 },
        .{ .star, 3, 3 },
        .{ .slash, 3, 3 },
        .{ .caret, 6, 5 },
    },
    .max_depth = 32,
});

const Src = struct {
    lx: Lexer,
    src: []const u8,
    cur: ?Tok = null,

    fn init(src: []const u8) Src {
        return .{ .lx = Lexer.init(src), .src = src };
    }
    pub fn peek(s: *Src) ?Tag {
        if (s.cur == null) s.cur = s.lx.next();
        return if (s.cur.?.tag == .eof) null else s.cur.?.tag;
    }
    pub fn next(s: *Src) Tok {
        const t = s.cur orelse s.lx.next();
        s.cur = null;
        return t;
    }
};

const Eval = struct {
    src: []const u8,
    pub const Node = f64;
    pub const Token = Tok;
    pub const Error = error{BadNumber};

    pub fn leaf(e: *Eval, t: Tok) Error!f64 {
        if (t.tag == .lparen) return error.BadNumber; // grouping not implemented
        return std.fmt.parseFloat(f64, t.text(e.src)) catch error.BadNumber;
    }
    pub fn unary(_: *Eval, op: Tag, x: f64) Error!f64 {
        return switch (op) {
            .minus => -x,
            else => unreachable,
        };
    }
    pub fn binary(_: *Eval, op: Tag, l: f64, r: f64) Error!f64 {
        return switch (op) {
            .plus => l + r,
            .minus => l - r,
            .star => l * r,
            .slash => l / r,
            .caret => std.math.pow(f64, l, r),
            else => unreachable,
        };
    }
};

fn evalStr(s: []const u8) !f64 {
    var src = Src.init(s);
    var e: Eval = .{ .src = s };
    return Parser.parse(&src, &e);
}

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    const cases = .{
        .{ "1 + 2 * 3 ^ 2", 19.0 },
        .{ "2 ^ 3 ^ 2", 512.0 },
        .{ "10 - 3 - 2", 5.0 },
        .{ "-3 * -4", 12.0 },
    };
    inline for (cases) |c| {
        const got = try evalStr(c[0]);
        try w.print("{s:>16} = {d:>7.1}\n", .{ c[0], got });
        if (got != c[1]) return error.ExampleFailed;
    }

    if (evalStr("1 +")) |_| return error.ExampleFailed else |err| {
        try w.print("{s:>16} -> {s}\n", .{ "1 +", @errorName(err) });
        if (err != error.UnexpectedEnd) return error.ExampleFailed;
    }

    // 31 unary minuses exceeds max_depth (32) before the leaf is even reached.
    var deep_buf: [64]u8 = undefined;
    for (&deep_buf) |*c| c.* = '-';
    deep_buf[63] = '1';
    if (evalStr(&deep_buf)) |_| return error.ExampleFailed else |err| {
        try w.print("{s:>16} -> {s}\n", .{ "63 unary minuses", @errorName(err) });
        if (err != error.TooDeep) return error.ExampleFailed;
    }
    try w.flush();
}
