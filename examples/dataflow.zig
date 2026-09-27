//! `graph.node.Graph`: build Const/Add/Accumulator nodes out of insertion
//! order, evaluate a few ticks, and show that a would-be cycle is rejected.
const std = @import("std");
const ts = @import("tsunami");
const node = ts.graph.node;

const Const = struct {
    pub const title = "const";
    pub const inputs = 0;
    pub const outputs = 1;
    v: f32,
    pub fn process(self: *Const, _: *const [0]f32, out: *[1]f32) void {
        out[0] = self.v;
    }
};

const Add = struct {
    pub const inputs = 2;
    pub const outputs = 1;
    pub fn process(_: *Add, in: *const [2]f32, out: *[1]f32) void {
        out[0] = in[0] + in[1];
    }
};

const Accumulator = struct {
    pub const inputs = 1;
    pub const outputs = 1;
    sum: f32 = 0,
    pub fn process(self: *Accumulator, in: *const [1]f32, out: *[1]f32) void {
        self.sum += in[0];
        out[0] = self.sum;
    }
};

const G = node.Graph(f32, .{ Const, Add, Accumulator });

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    var gpa_state: std.heap.DebugAllocator(.{}) = .init;
    defer std.debug.assert(gpa_state.deinit() == .ok);
    const gpa = gpa_state.allocator();

    var g: G = .empty;
    defer g.deinit(gpa);

    // Inserted out of dependency order: the accumulator and adder are added
    // before their sources, but `evaluate` still runs in topological order.
    const acc = try g.add(gpa, Accumulator{});
    const add = try g.add(gpa, Add{});
    const a = try g.add(gpa, Const{ .v = 2 });
    const b = try g.add(gpa, Const{ .v = 3 });
    try g.connect(gpa, .{ .node = a, .port = 0 }, .{ .node = add, .port = 0 });
    try g.connect(gpa, .{ .node = b, .port = 0 }, .{ .node = add, .port = 1 });
    try g.connect(gpa, .{ .node = add, .port = 0 }, .{ .node = acc, .port = 0 });

    var running: f32 = 0;
    for (0..4) |tick| {
        try g.evaluate(gpa);
        const out = try g.output(acc, 0);
        running += 5;
        try w.print("tick {d}: accumulator = {d:.0} (want {d:.0})\n", .{ tick, out, running });
        if (out != running) return error.ExampleFailed;
    }

    // add -> acc -> add would close a cycle through the already-connected edge.
    if (g.connect(gpa, .{ .node = acc, .port = 0 }, .{ .node = add, .port = 0 })) |_| {
        return error.ExampleFailed;
    } else |err| {
        try w.print("cycle rejected: {s}\n", .{@errorName(err)});
        if (err != error.Cycle) return error.ExampleFailed;
    }
    try w.flush();
}
