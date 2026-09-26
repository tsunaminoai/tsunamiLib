const std = @import("std");

pub const Error = error{ UnknownFlag, MissingValue, BadValue, DuplicateFlag };

/// Where parsing failed, for error messages. Index into argv.
pub const Diag = struct { index: usize = 0, arg: []const u8 = "" };

pub fn Result(comptime Opts: type) type {
    return struct { opts: Opts, rest: []const []const u8 };
}

/// Parse argv (without argv[0]) into `Opts`, a struct whose fields are the
/// flags: `foo_bar` ↔ `--foo-bar`. Field types: bool (switch; `--no-x`
/// clears), ints (base prefixes allowed), floats, enums, `[]const u8`, and
/// `?T` of those. Fields without defaults are required. Optional decls:
/// `pub const short = .{ .v = .verbose };` and `pub const help = .{ .verbose = "..." };`.
/// Parsing stops at `--` or the first positional; the tail is `rest`.
/// Unknown flags are errors — a typo must never fall back to a default.
pub fn parse(comptime Opts: type, argv: []const []const u8, diag: ?*Diag) Error!Result(Opts) {
    const fields = @typeInfo(Opts).@"struct".fields;
    var out: Opts = undefined;
    var seen: std.StaticBitSet(fields.len) = .initEmpty();

    var i: usize = 0;
    while (i < argv.len) : (i += 1) {
        const arg = argv[i];
        if (std.mem.eql(u8, arg, "--")) {
            i += 1;
            break;
        }
        if (arg.len < 2 or arg[0] != '-') break;
        errdefer if (diag) |d| {
            d.* = .{ .index = i, .arg = arg };
        };

        var name: []const u8 = undefined;
        var inline_val: ?[]const u8 = null;
        if (arg[1] == '-') {
            const body = arg[2..];
            if (std.mem.indexOfScalar(u8, body, '=')) |eq| {
                name = body[0..eq];
                inline_val = body[eq + 1 ..];
            } else name = body;
        } else {
            if (arg.len != 2) return error.UnknownFlag;
            name = shortName(Opts, arg[1]) orelse return error.UnknownFlag;
        }

        var negate = false;
        const idx = lookup(Opts).get(name) orelse blk: {
            if (std.mem.startsWith(u8, name, "no-")) if (lookup(Opts).get(name[3..])) |j| {
                negate = true;
                break :blk j;
            };
            return error.UnknownFlag;
        };

        inline for (fields, 0..) |f, fi| if (fi == idx) {
            if (seen.isSet(fi)) return error.DuplicateFlag;
            seen.set(fi);
            if (comptime isBool(f.type)) {
                if (inline_val != null) return error.BadValue;
                @field(out, f.name) = !negate;
            } else {
                if (negate) return error.UnknownFlag;
                const v = inline_val orelse v: {
                    i += 1;
                    if (i >= argv.len) return error.MissingValue;
                    break :v argv[i];
                };
                @field(out, f.name) = try parseValue(f.type, v);
            }
        };
    }

    inline for (fields, 0..) |f, fi| if (!seen.isSet(fi)) {
        if (f.defaultValue()) |d| {
            @field(out, f.name) = d;
        } else if (@typeInfo(f.type) == .optional) {
            @field(out, f.name) = null;
        } else if (comptime isBool(f.type)) {
            @field(out, f.name) = false;
        } else {
            if (diag) |d| d.* = .{ .index = argv.len, .arg = comptime flagName(f.name) };
            return error.MissingValue;
        }
    };
    return .{ .opts = out, .rest = argv[i..] };
}

/// Comptime usage text: one line per flag with type, default and help.
pub fn usage(comptime Opts: type) []const u8 {
    comptime {
        var s: []const u8 = "";
        for (@typeInfo(Opts).@"struct".fields) |f| {
            var line: []const u8 = "  ";
            if (@hasDecl(Opts, "short")) for (@typeInfo(@TypeOf(Opts.short)).@"struct".fields) |sf| {
                if (std.mem.eql(u8, @tagName(@field(Opts.short, sf.name)), f.name)) line = line ++ "-" ++ sf.name ++ ", ";
            };
            line = line ++ "--" ++ flagName(f.name);
            if (!isBool(f.type)) line = line ++ " <" ++ typeLabel(f.type) ++ ">";
            if (@hasDecl(Opts, "help") and @hasField(@TypeOf(Opts.help), f.name)) line = line ++ "  " ++ @field(Opts.help, f.name);
            if (f.defaultValue()) |d| if (!isBool(f.type)) {
                line = line ++ " (default: " ++ fmtDefault(f.type, d) ++ ")";
            };
            if (f.default_value_ptr == null and @typeInfo(f.type) != .optional and !isBool(f.type)) line = line ++ " (required)";
            s = s ++ line ++ "\n";
        }
        return s;
    }
}

// ── Internals ────────────────────────────────────────────────────────────

fn lookup(comptime Opts: type) type {
    const fields = @typeInfo(Opts).@"struct".fields;
    const KV = struct { []const u8, usize };
    var kvs: [fields.len]KV = undefined;
    for (fields, 0..) |f, i| kvs[i] = .{ flagName(f.name), i };
    const final = kvs;
    return struct {
        const map = std.StaticStringMap(usize).initComptime(final);
        fn get(name: []const u8) ?usize {
            return map.get(name);
        }
    };
}

fn shortName(comptime Opts: type, c: u8) ?[]const u8 {
    if (!@hasDecl(Opts, "short")) return null;
    inline for (@typeInfo(@TypeOf(Opts.short)).@"struct".fields) |sf| {
        comptime std.debug.assert(sf.name.len == 1);
        if (sf.name[0] == c) return comptime flagName(@tagName(@field(Opts.short, sf.name)));
    }
    return null;
}

fn flagName(comptime field: []const u8) []const u8 {
    comptime {
        var buf: [field.len]u8 = undefined;
        for (field, 0..) |ch, i| buf[i] = if (ch == '_') '-' else ch;
        const final = buf;
        return &final;
    }
}

fn isBool(comptime T: type) bool {
    return T == bool or T == ?bool;
}

fn parseValue(comptime T: type, v: []const u8) Error!T {
    return switch (@typeInfo(T)) {
        .optional => |o| try parseValue(o.child, v),
        .int => std.fmt.parseInt(T, v, 0) catch error.BadValue,
        .float => std.fmt.parseFloat(T, v) catch error.BadValue,
        .@"enum" => std.meta.stringToEnum(T, v) orelse error.BadValue,
        .pointer => if (T == []const u8) v else @compileError("unsupported flag type " ++ @typeName(T)),
        else => @compileError("unsupported flag type " ++ @typeName(T)),
    };
}

fn typeLabel(comptime T: type) []const u8 {
    return switch (@typeInfo(T)) {
        .optional => |o| typeLabel(o.child),
        .int => "int",
        .float => "num",
        .@"enum" => |e| blk: {
            var s: []const u8 = "";
            for (e.fields, 0..) |ef, i| s = s ++ (if (i == 0) "" else "|") ++ ef.name;
            break :blk s;
        },
        else => "str",
    };
}

fn fmtDefault(comptime T: type, comptime d: T) []const u8 {
    return switch (@typeInfo(T)) {
        .@"enum" => @tagName(d),
        .pointer => d,
        .optional => if (d) |x| fmtDefault(@typeInfo(T).optional.child, x) else "none",
        else => std.fmt.comptimePrint("{d}", .{d}),
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

const TOpts = struct {
    matches: u32 = 96,
    floor: []const u8 = "default",
    fast: bool = false,
    ratio: f32 = 0.5,
    mode: enum { a, b } = .a,
    seed: ?u64 = null,
    out_dir: []const u8,

    pub const short = .{ .f = .fast, .o = .out_dir };
    pub const help = .{ .matches = "games to play" };
};

test "parses flags, defaults, rest" {
    const argv = [_][]const u8{ "--matches", "0x10", "-f", "--out-dir=/tmp", "--mode", "b", "file1", "--fast" };
    const r = try parse(TOpts, &argv, null);
    try testing.expectEqual(@as(u32, 16), r.opts.matches);
    try testing.expect(r.opts.fast);
    try testing.expectEqualStrings("/tmp", r.opts.out_dir);
    try testing.expectEqualStrings("default", r.opts.floor);
    try testing.expectEqual(.b, r.opts.mode);
    try testing.expectEqual(@as(?u64, null), r.opts.seed);
    try testing.expectEqual(@as(usize, 2), r.rest.len);
}

test "errors carry diagnostics" {
    var d: Diag = .{};
    try testing.expectError(error.UnknownFlag, parse(TOpts, &.{ "-o", "x", "--matchez", "1" }, &d));
    try testing.expectEqual(@as(usize, 2), d.index);
    try testing.expectError(error.MissingValue, parse(TOpts, &.{ "-o", "x", "--matches" }, &d));
    try testing.expectError(error.BadValue, parse(TOpts, &.{ "-o", "x", "--matches", "lots" }, &d));
    try testing.expectError(error.BadValue, parse(TOpts, &.{ "-o", "x", "--mode", "z" }, &d));
    try testing.expectError(error.DuplicateFlag, parse(TOpts, &.{ "-o", "x", "-o", "y" }, &d));
    try testing.expectError(error.MissingValue, parse(TOpts, &.{"--fast"}, &d));
    try testing.expectEqualStrings("out-dir", d.arg);
}

test "--no- negation and -- terminator" {
    const r = try parse(TOpts, &.{ "-o", "x", "--no-fast", "--", "--matches" }, null);
    try testing.expect(!r.opts.fast);
    try testing.expectEqualStrings("--matches", r.rest[0]);
}

test "usage is comptime" {
    const u = comptime usage(TOpts);
    try testing.expect(std.mem.indexOf(u8, u, "--matches <int>  games to play (default: 96)") != null);
    try testing.expect(std.mem.indexOf(u8, u, "-o, --out-dir <str> (required)") != null);
    try testing.expect(std.mem.indexOf(u8, u, "--mode <a|b>") != null);
}

test "fuzz parse" {
    try testing.fuzz({}, struct {
        fn f(_: void, s: *testing.Smith) !void {
            var buf: [256]u8 = undefined;
            const n = s.slice(&buf);
            var parts: [16][]const u8 = undefined;
            var np: usize = 0;
            var it = std.mem.splitScalar(u8, buf[0..n], 0);
            while (it.next()) |p| : (np += 1) {
                if (np == parts.len) break;
                parts[np] = p;
            }
            _ = parse(TOpts, parts[0..np], null) catch return;
        }
    }.f, .{});
}
