//! text.args.parse over a flags struct: short flags, enum, optional,
//! required field, comptime usage text, and the Diag error path.
//! With no real argv (as under `zig build examples`) it parses a built-in
//! demo argv instead, so the example is always runnable and deterministic.
const std = @import("std");
const ts = @import("tsunami");
const args = ts.text.args;

const Opts = struct {
    port: u16 = 8080,
    mode: enum { serve, bench } = .serve,
    tag: ?[]const u8 = null,
    verbose: bool = false,
    name: []const u8,

    pub const short = .{ .v = .verbose, .p = .port };
    pub const help = .{
        .port = "listen port",
        .mode = "operating mode",
        .name = "instance name",
    };
};

const demo_argv = [_][]const u8{ "-p", "9090", "--mode", "bench", "--name", "spectrum-node", "-v" };
const bad_argv = [_][]const u8{ "--name", "x", "--modee", "bench" };

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    const real_argv = try init.minimal.args.toSlice(init.arena.allocator());
    const argv: []const []const u8 = if (real_argv.len > 1) real_argv[1..] else &demo_argv;

    const result = try args.parse(Opts, argv, null);
    try w.print("parsed: name={s} port={d} mode={s} tag={?s} verbose={}\n", .{
        result.opts.name, result.opts.port, @tagName(result.opts.mode), result.opts.tag, result.opts.verbose,
    });

    try w.print("\nusage:\n{s}", .{comptime args.usage(Opts)});

    var diag: args.Diag = .{};
    const bad = args.parse(Opts, &bad_argv, &diag);
    try w.print("\nbad argv --modee: {s} at argv[{d}] (\"{s}\")\n", .{
        if (bad) |_| "accepted (BUG)" else |e| @errorName(e),
        diag.index,
        diag.arg,
    });
    try w.flush();

    const ok = std.mem.eql(u8, result.opts.name, "spectrum-node") and
        result.opts.port == 9090 and result.opts.mode == .bench and result.opts.verbose and
        bad == error.UnknownFlag;
    if (!ok) return error.ExampleFailed;
}
