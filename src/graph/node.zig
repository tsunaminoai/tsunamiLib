const std = @import("std");

pub const Error = error{ BadNode, BadPort, Cycle } || std.mem.Allocator.Error;

pub const Id = enum(u32) { _ };

/// Dataflow graph over value type `V`. Each kind in `kinds` is a struct whose
/// fields are its state, declaring `pub const inputs`, `pub const outputs`
/// and `pub fn process(*Self, *const [inputs]V, *[outputs]V) void`.
///
/// Nodes live in a comptime-built tagged union and evaluation dispatches with
/// `inline else`, so every `process` call is static and inlinable — no
/// vtables, no type-erased state blobs.
pub fn Graph(comptime V: type, comptime kinds: anytype) type {
    const n_kinds = kinds.len;
    comptime var names: [n_kinds][]const u8 = undefined;
    comptime var types: [n_kinds]type = undefined;
    comptime var max_in: usize = 0;
    comptime var max_out: usize = 0;
    inline for (kinds, 0..) |K, i| {
        if (!@hasDecl(K, "inputs") or !@hasDecl(K, "outputs") or !@hasDecl(K, "process"))
            @compileError(@typeName(K) ++ " needs inputs, outputs, process");
        names[i] = if (@hasDecl(K, "title")) K.title else std.fmt.comptimePrint("k{d}", .{i});
        types[i] = K;
        max_in = @max(max_in, K.inputs);
        max_out = @max(max_out, K.outputs);
    }
    const KindT = @Enum(u16, .exhaustive, &names, &std.simd.iota(u16, n_kinds));
    const NodeT = @Union(.auto, KindT, &names, &types, &@splat(.{}));

    return struct {
        const Self = @This();
        pub const Kind = KindT;
        pub const Node = NodeT;
        pub const Port = struct { node: Id, port: u8 };

        const Slot = struct {
            node: Node,
            alive: bool = true,
            src: [max_in]?Port = @splat(null),
            out: [max_out]V = @splat(std.mem.zeroes(V)),
        };

        slots: std.ArrayList(Slot) = .empty,
        order: std.ArrayList(u32) = .empty,
        dirty: bool = false,

        pub const empty: Self = .{};

        pub fn deinit(g: *Self, gpa: std.mem.Allocator) void {
            g.slots.deinit(gpa);
            g.order.deinit(gpa);
            g.* = undefined;
        }

        pub fn add(g: *Self, gpa: std.mem.Allocator, node: anytype) Error!Id {
            const T = @TypeOf(node);
            const tag = comptime kindOf(T);
            try g.slots.append(gpa, .{ .node = @unionInit(Node, @tagName(tag), node) });
            g.dirty = true;
            return @enumFromInt(g.slots.items.len - 1);
        }

        pub fn remove(g: *Self, id: Id) Error!void {
            const s = try g.slot(id);
            s.alive = false;
            for (g.slots.items) |*o| for (&o.src) |*p| {
                if (p.* != null and p.*.?.node == id) p.* = null;
            };
            g.dirty = true;
        }

        pub fn get(g: *Self, comptime T: type, id: Id) Error!*T {
            const s = try g.slot(id);
            return switch (s.node) {
                comptime kindOf(T) => |*n| n,
                else => error.BadNode,
            };
        }

        pub fn output(g: *const Self, id: Id, port: u8) Error!V {
            const i = @intFromEnum(id);
            if (i >= g.slots.items.len or !g.slots.items[i].alive) return error.BadNode;
            const s = &g.slots.items[i];
            if (port >= portCount(s.node, .out)) return error.BadPort;
            return s.out[port];
        }

        /// Rejects edges that would close a cycle; a replaced input is overwritten.
        pub fn connect(g: *Self, gpa: std.mem.Allocator, from: Port, to: Port) Error!void {
            const src = try g.slot(from.node);
            const dst = try g.slot(to.node);
            if (from.port >= portCount(src.node, .out) or to.port >= portCount(dst.node, .in)) return error.BadPort;
            if (from.node == to.node or try g.reaches(gpa, to.node, from.node)) return error.Cycle;
            dst.src[to.port] = from;
            g.dirty = true;
        }

        pub fn disconnect(g: *Self, to: Port) Error!void {
            const dst = try g.slot(to.node);
            if (to.port >= portCount(dst.node, .in)) return error.BadPort;
            dst.src[to.port] = null;
            g.dirty = true;
        }

        /// One pass in topological order. Unconnected inputs read zero.
        pub fn evaluate(g: *Self, gpa: std.mem.Allocator) Error!void {
            if (g.dirty) try g.sort(gpa);
            const items = g.slots.items;
            for (g.order.items) |i| {
                const s = &items[i];
                switch (s.node) {
                    inline else => |*n| {
                        const K = @TypeOf(n.*);
                        var in: [K.inputs]V = undefined;
                        inline for (&in, 0..) |*v, p| v.* = if (s.src[p]) |sp| items[@intFromEnum(sp.node)].out[sp.port] else std.mem.zeroes(V);
                        n.process(&in, s.out[0..K.outputs]);
                    },
                }
            }
        }

        // ── Internals ────────────────────────────────────────────────────

        fn kindOf(comptime T: type) Kind {
            inline for (kinds, 0..) |K, i| if (K == T) return @enumFromInt(i);
            @compileError(@typeName(T) ++ " is not a registered node kind");
        }

        fn portCount(n: Node, comptime dir: enum { in, out }) usize {
            return switch (n) {
                inline else => |_, tag| blk: {
                    const K = @FieldType(Node, @tagName(tag));
                    break :blk if (dir == .in) K.inputs else K.outputs;
                },
            };
        }

        fn slot(g: *Self, id: Id) Error!*Slot {
            const i = @intFromEnum(id);
            if (i >= g.slots.items.len or !g.slots.items[i].alive) return error.BadNode;
            return &g.slots.items[i];
        }

        /// Does `start` feed (transitively) into `target`? Walks upstream from target.
        fn reaches(g: *const Self, gpa: std.mem.Allocator, start: Id, target: Id) Error!bool {
            var visited = try std.DynamicBitSetUnmanaged.initEmpty(gpa, g.slots.items.len);
            defer visited.deinit(gpa);
            var stack: std.ArrayList(Id) = .empty;
            defer stack.deinit(gpa);
            try stack.append(gpa, target);
            while (stack.pop()) |cur| {
                if (cur == start) return true;
                const ci = @intFromEnum(cur);
                if (visited.isSet(ci)) continue;
                visited.set(ci);
                for (g.slots.items[ci].src) |p| if (p) |pp| try stack.append(gpa, pp.node);
            }
            return false;
        }

        fn sort(g: *Self, gpa: std.mem.Allocator) Error!void {
            const items = g.slots.items;
            g.order.clearRetainingCapacity();
            try g.order.ensureTotalCapacity(gpa, items.len);
            const done = try gpa.alloc(bool, items.len);
            defer gpa.free(done);
            @memset(done, false);
            // Repeated relaxation: acyclicity is enforced by connect(), so each
            // sweep places at least one node. O(n²) worst case, run only when dirty.
            var placed: usize = 0;
            var live: usize = 0;
            for (items) |s| live += @intFromBool(s.alive);
            while (placed < live) {
                const before = placed;
                for (items, 0..) |s, i| {
                    if (!s.alive or done[i]) continue;
                    const ready = for (s.src) |p| {
                        if (p) |pp| if (!done[@intFromEnum(pp.node)]) break false;
                    } else true;
                    if (ready) {
                        done[i] = true;
                        g.order.appendAssumeCapacity(@intCast(i));
                        placed += 1;
                    }
                }
                std.debug.assert(placed > before);
            }
            g.dirty = false;
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

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

const Accum = struct {
    pub const inputs = 1;
    pub const outputs = 1;
    sum: f32 = 0,
    pub fn process(self: *Accum, in: *const [1]f32, out: *[1]f32) void {
        self.sum += in[0];
        out[0] = self.sum;
    }
};

const G = Graph(f32, .{ Const, Add, Accum });

test Graph {
    var g: G = .empty;
    defer g.deinit(testing.allocator);
    const acc = try g.add(testing.allocator, Accum{});
    const add = try g.add(testing.allocator, Add{});
    const a = try g.add(testing.allocator, Const{ .v = 2 });
    const b = try g.add(testing.allocator, Const{ .v = 3 });
    try g.connect(testing.allocator, .{ .node = a, .port = 0 }, .{ .node = add, .port = 0 });
    try g.connect(testing.allocator, .{ .node = b, .port = 0 }, .{ .node = add, .port = 1 });
    try g.connect(testing.allocator, .{ .node = add, .port = 0 }, .{ .node = acc, .port = 0 });
    try g.evaluate(testing.allocator);
    try g.evaluate(testing.allocator);
    try testing.expectEqual(@as(f32, 10), try g.output(acc, 0));
    (try g.get(Const, a)).v = 10;
    try g.evaluate(testing.allocator);
    try testing.expectEqual(@as(f32, 23), try g.output(acc, 0));
}

test "cycles, bad ports and removed nodes are rejected" {
    var g: G = .empty;
    defer g.deinit(testing.allocator);
    const x = try g.add(testing.allocator, Accum{});
    const y = try g.add(testing.allocator, Accum{});
    try g.connect(testing.allocator, .{ .node = x, .port = 0 }, .{ .node = y, .port = 0 });
    try testing.expectError(error.Cycle, g.connect(testing.allocator, .{ .node = y, .port = 0 }, .{ .node = x, .port = 0 }));
    try testing.expectError(error.Cycle, g.connect(testing.allocator, .{ .node = x, .port = 0 }, .{ .node = x, .port = 0 }));
    try testing.expectError(error.BadPort, g.connect(testing.allocator, .{ .node = x, .port = 1 }, .{ .node = y, .port = 0 }));
    try testing.expectError(error.BadNode, g.get(Const, x));
    try g.remove(x);
    try testing.expectError(error.BadNode, g.output(x, 0));
    try g.evaluate(testing.allocator);
    try testing.expectEqual(@as(f32, 0), try g.output(y, 0));
}
