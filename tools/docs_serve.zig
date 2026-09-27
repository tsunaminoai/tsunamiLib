//! Serves zig-out/docs on 127.0.0.1 (autodoc can't load over file://).
//! Usage: docs_serve <docs-dir> [port]. Only the four autodoc files are
//! served — no path handling, so no traversal surface.
const std = @import("std");
const Io = std.Io;

const files = [_]struct { path: []const u8, name: []const u8, mime: []const u8 }{
    .{ .path = "/", .name = "index.html", .mime = "text/html" },
    .{ .path = "/main.js", .name = "main.js", .mime = "application/javascript" },
    .{ .path = "/main.wasm", .name = "main.wasm", .mime = "application/wasm" },
    .{ .path = "/sources.tar", .name = "sources.tar", .mime = "application/x-tar" },
};

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const arena = init.arena.allocator();
    const argv = try init.minimal.args.toSlice(arena);
    if (argv.len < 2) return error.Usage;
    const port = if (argv.len > 2) try std.fmt.parseInt(u16, argv[2], 10) else 8080;

    var dir = try Io.Dir.cwd().openDir(io, argv[1], .{});
    defer dir.close(io);
    var bodies: [files.len][]const u8 = undefined;
    for (files, &bodies) |f, *b| b.* = try dir.readFileAlloc(io, f.name, arena, .limited(256 << 20));

    const address = try Io.net.IpAddress.parse("127.0.0.1", port);
    var server = try address.listen(io, .{ .reuse_address = true });
    defer server.deinit(io);
    std.log.info("docs at http://127.0.0.1:{d}/ (ctrl-c to stop)", .{server.socket.address.getPort()});

    while (true) {
        const stream = server.accept(io) catch continue;
        defer stream.close(io);
        serve(io, stream, &bodies) catch |e| std.log.warn("connection: {t}", .{e});
    }
}

fn serve(io: Io, stream: Io.net.Stream, bodies: *const [files.len][]const u8) !void {
    var rbuf: [8192]u8 = undefined;
    var wbuf: [8192]u8 = undefined;
    var r = stream.reader(io, &rbuf);
    var w = stream.writer(io, &wbuf);
    var http: std.http.Server = .init(&r.interface, &w.interface);
    while (true) {
        var req = http.receiveHead() catch |e| switch (e) {
            error.HttpConnectionClosing => return,
            else => return e,
        };
        const t = req.head.target;
        const path = t[0 .. std.mem.indexOfAny(u8, t, "?#") orelse t.len];
        const i = for (files, 0..) |f, i| {
            if (std.mem.eql(u8, f.path, path)) break i;
        } else if (std.mem.eql(u8, path, "/index.html")) 0 else {
            try req.respond("not found", .{ .status = .not_found });
            continue;
        };
        try req.respond(bodies[i], .{ .extra_headers = &.{
            .{ .name = "content-type", .value = files[i].mime },
            .{ .name = "cache-control", .value = "no-cache" },
        } });
    }
}
