const std = @import("std");

/// Binding powers (≥ 1): higher binds tighter. Infix `left == right` is
/// left-associative; `left > right` is right-associative.
pub fn Spec(comptime Tag: type) type {
    return struct {
        prefix: []const struct { Tag, u8 } = &.{},
        infix: []const struct { Tag, u8, u8 } = &.{},
        postfix: []const struct { Tag, u8 } = &.{},
        max_depth: u16 = 256,
    };
}

/// Precedence-climbing driver. Tables are flattened at comptime into arrays
/// indexed by tag, so each lookup is a constant-offset load.
///
/// `src` needs `peek() ?Tag` (null at end) and `next() Token`.
/// `b` needs `Node`, `Token`, `Error`, and `leaf(Token) Error!Node`,
/// `unary(Tag, Node) Error!Node`, `binary(Tag, Node, Node) Error!Node`.
/// Recursion is capped at `max_depth`, so hostile input can't blow the stack.
pub fn Parser(comptime Tag: type, comptime spec: Spec(Tag)) type {
    const n = @typeInfo(Tag).@"enum".fields.len;
    const NoBp = std.math.maxInt(u8);
    const Tbl = struct { pre: [n]u8, in_l: [n]u8, in_r: [n]u8, post: [n]u8 };
    const tbl: Tbl = blk: {
        var t: Tbl = .{ .pre = @splat(NoBp), .in_l = @splat(NoBp), .in_r = @splat(NoBp), .post = @splat(NoBp) };
        for (spec.prefix) |e| t.pre[@intFromEnum(e[0])] = e[1];
        for (spec.infix) |e| {
            t.in_l[@intFromEnum(e[0])] = e[1];
            t.in_r[@intFromEnum(e[0])] = e[2];
        }
        for (spec.postfix) |e| t.post[@intFromEnum(e[0])] = e[1];
        break :blk t;
    };

    return struct {
        pub fn Error(comptime B: type) type {
            return error{ UnexpectedEnd, TooDeep } || B.Error;
        }

        pub fn parse(src: anytype, b: anytype) Error(@TypeOf(b.*))!@TypeOf(b.*).Node {
            return expr(src, b, 0, 0);
        }

        fn expr(src: anytype, b: anytype, min_bp: u8, depth: u16) Error(@TypeOf(b.*))!@TypeOf(b.*).Node {
            if (depth >= spec.max_depth) return error.TooDeep;
            const first = src.peek() orelse return error.UnexpectedEnd;
            var lhs = if (tbl.pre[@intFromEnum(first)] != NoBp) blk: {
                _ = src.next();
                const x = try expr(src, b, tbl.pre[@intFromEnum(first)], depth + 1);
                break :blk try b.unary(first, x);
            } else try b.leaf(src.next());

            while (src.peek()) |tag| {
                const i = @intFromEnum(tag);
                if (tbl.post[i] != NoBp) {
                    if (tbl.post[i] <= min_bp) break;
                    _ = src.next();
                    lhs = try b.unary(tag, lhs);
                } else if (tbl.in_l[i] != NoBp) {
                    if (tbl.in_l[i] <= min_bp) break;
                    _ = src.next();
                    const rhs = try expr(src, b, tbl.in_r[i], depth + 1);
                    lhs = try b.binary(tag, lhs, rhs);
                } else break;
            }
            return lhs;
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

const ATag = enum { num, plus, minus, star, caret, bang };
const Tok = struct { tag: ATag, v: f64 = 0 };

const Src = struct {
    toks: []const Tok,
    i: usize = 0,
    fn peek(s: *Src) ?ATag {
        return if (s.i < s.toks.len) s.toks[s.i].tag else null;
    }
    fn next(s: *Src) Tok {
        defer s.i += 1;
        return s.toks[s.i];
    }
};

const Eval = struct {
    pub const Node = f64;
    pub const Token = Tok;
    pub const Error = error{NotANumber};
    fn leaf(_: *Eval, t: Tok) Error!f64 {
        return if (t.tag == .num) t.v else error.NotANumber;
    }
    fn unary(_: *Eval, op: ATag, x: f64) Error!f64 {
        return switch (op) {
            .minus => -x,
            .bang => std.math.gamma(f64, x + 1),
            else => unreachable,
        };
    }
    fn binary(_: *Eval, op: ATag, l: f64, r: f64) Error!f64 {
        return switch (op) {
            .plus => l + r,
            .minus => l - r,
            .star => l * r,
            .caret => std.math.pow(f64, l, r),
            else => unreachable,
        };
    }
};

const P = Parser(ATag, .{
    .prefix = &.{.{ .minus, 7 }},
    .infix = &.{ .{ .plus, 1, 1 }, .{ .minus, 1, 1 }, .{ .star, 3, 3 }, .{ .caret, 6, 5 } },
    .postfix = &.{.{ .bang, 9 }},
    .max_depth = 32,
});

fn run(toks: []const Tok) !f64 {
    var s: Src = .{ .toks = toks };
    var e: Eval = .{};
    return P.parse(&s, &e);
}

fn num(v: f64) Tok {
    return .{ .tag = .num, .v = v };
}

test "precedence and associativity" {
    // 1 + 2 * 3 ^ 2 = 19
    try testing.expectEqual(@as(f64, 19), try run(&.{ num(1), .{ .tag = .plus }, num(2), .{ .tag = .star }, num(3), .{ .tag = .caret }, num(2) }));
    // 2 ^ 3 ^ 2 = 512 (right-assoc)
    try testing.expectEqual(@as(f64, 512), try run(&.{ num(2), .{ .tag = .caret }, num(3), .{ .tag = .caret }, num(2) }));
    // 10 - 3 - 2 = 5 (left-assoc)
    try testing.expectEqual(@as(f64, 5), try run(&.{ num(10), .{ .tag = .minus }, num(3), .{ .tag = .minus }, num(2) }));
    // -3! = -6 (postfix binds tighter than prefix)
    try testing.expectApproxEqAbs(@as(f64, -6), try run(&.{ .{ .tag = .minus }, num(3), .{ .tag = .bang } }), 1e-9);
}

test "errors: dangling operator, bad leaf, depth cap" {
    try testing.expectError(error.UnexpectedEnd, run(&.{ num(1), .{ .tag = .plus } }));
    try testing.expectError(error.UnexpectedEnd, run(&.{}));
    try testing.expectError(error.NotANumber, run(&.{.{ .tag = .star }}));
    var deep: [64]Tok = @splat(.{ .tag = .minus });
    deep[63] = num(1);
    try testing.expectError(error.TooDeep, run(&deep));
}

test "fuzz parse never crashes" {
    try testing.fuzz({}, struct {
        fn f(_: void, s: *testing.Smith) !void {
            var buf: [64]u8 = undefined;
            const len = s.slice(&buf);
            var toks: [64]Tok = undefined;
            for (buf[0..len], toks[0..len]) |c, *t| t.* = .{ .tag = @enumFromInt(c % 6), .v = @floatFromInt(c) };
            _ = run(toks[0..len]) catch return;
        }
    }.f, .{});
}
