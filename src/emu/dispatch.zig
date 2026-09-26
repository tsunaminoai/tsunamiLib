const std = @import("std");

pub const Stop = enum { cont, halt };

pub const Error = error{IllegalOpcode};

/// Comptime-built 8-bit opcode ISA. `entries` is a tuple of
/// `.{ .code = 0x69, .op = .adc, .mode = .imm, .cycles = 2 }`; `.op` names a
/// decl on `Cpu` with signature `fn (*Cpu, comptime Mode) void|Stop`.
///
/// `step`/`run` switch with `inline else`, so each opcode prong calls its
/// handler with a comptime-known mode: addressing is monomorphised per opcode
/// and `run` compiles to direct-threaded code via labeled `continue :sw`.
pub fn Isa(comptime Cpu: type, comptime Mode: type, comptime entries: anytype) type {
    const Spec = struct { op: @EnumLiteral(), mode: Mode, cycles: u8 };
    const specs: [256]?Spec = blk: {
        var t: [256]?Spec = @splat(null);
        for (entries) |e| {
            if (t[e.code] != null) @compileError(std.fmt.comptimePrint("duplicate opcode 0x{X:0>2}", .{e.code}));
            if (!@hasDecl(Cpu, @tagName(e.op))) @compileError("Cpu has no handler '" ++ @tagName(e.op) ++ "'");
            t[e.code] = .{ .op = e.op, .mode = e.mode, .cycles = e.cycles };
        }
        break :blk t;
    };

    return struct {
        pub const Info = struct { name: []const u8, mode: Mode, cycles: u8 };

        /// Runtime-queryable table for disassemblers/debuggers.
        pub const info: [256]?Info = blk: {
            var t: [256]?Info = @splat(null);
            for (specs, 0..) |s, i| if (s) |v| {
                t[i] = .{ .name = @tagName(v.op), .mode = v.mode, .cycles = v.cycles };
            };
            break :blk t;
        };

        inline fn exec(cpu: *Cpu, comptime code: u8) Stop {
            const s = comptime specs[code].?;
            const r = @field(Cpu, @tagName(s.op))(cpu, s.mode);
            return if (@TypeOf(r) == Stop) r else .cont;
        }

        /// Execute one opcode; returns its base cycle count.
        pub fn step(cpu: *Cpu, opcode: u8) Error!struct { cycles: u8, stop: Stop } {
            switch (opcode) {
                inline else => |c| {
                    if (comptime specs[c] == null) return error.IllegalOpcode;
                    return .{ .cycles = comptime specs[c].?.cycles, .stop = exec(cpu, c) };
                },
            }
        }

        /// Fetch–execute until a handler halts or `budget` cycles elapse.
        /// `Cpu.fetch(*Cpu) u8` supplies the next opcode. Returns cycles used.
        pub fn run(cpu: *Cpu, budget: u64) Error!u64 {
            var used: u64 = 0;
            sw: switch (cpu.fetch()) {
                inline else => |c| {
                    if (comptime specs[c] == null) return error.IllegalOpcode;
                    used += comptime specs[c].?.cycles;
                    if (exec(cpu, c) == .halt or used >= budget) return used;
                    continue :sw cpu.fetch();
                },
            }
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

const Toy = struct {
    const Mode = enum { imp, imm, zp };
    mem: [256]u8 = @splat(0),
    pc: u8 = 0,
    a: u8 = 0,

    fn fetch(self: *Toy) u8 {
        defer self.pc +%= 1;
        return self.mem[self.pc];
    }

    fn operand(self: *Toy, comptime mode: Mode) u8 {
        return switch (mode) {
            .imp => unreachable,
            .imm => self.fetch(),
            .zp => self.mem[self.fetch()],
        };
    }

    pub fn inc(self: *Toy, comptime _: Mode) void {
        self.a +%= 1;
    }
    pub fn add(self: *Toy, comptime mode: Mode) void {
        self.a +%= self.operand(mode);
    }
    pub fn hlt(_: *Toy, comptime _: Mode) Stop {
        return .halt;
    }

    const isa = Isa(Toy, Mode, .{
        .{ .code = 0x01, .op = .inc, .mode = .imp, .cycles = 1 },
        .{ .code = 0x02, .op = .add, .mode = .imm, .cycles = 2 },
        .{ .code = 0x03, .op = .add, .mode = .zp, .cycles = 3 },
        .{ .code = 0xFF, .op = .hlt, .mode = .imp, .cycles = 1 },
    });
};

test "run executes program until halt" {
    var cpu: Toy = .{};
    cpu.mem[0x80] = 10;
    @memcpy(cpu.mem[0..7], &[_]u8{ 0x01, 0x02, 5, 0x03, 0x80, 0x01, 0xFF });
    const cycles = try Toy.isa.run(&cpu, 1000);
    try testing.expectEqual(@as(u8, 17), cpu.a);
    try testing.expectEqual(@as(u64, 1 + 2 + 3 + 1 + 1), cycles);
}

test "illegal opcode and budget" {
    var cpu: Toy = .{};
    cpu.mem[0] = 0x42;
    try testing.expectError(error.IllegalOpcode, Toy.isa.run(&cpu, 10));
    cpu = .{};
    @memset(cpu.mem[0..], 0x01);
    try testing.expectEqual(@as(u64, 5), try Toy.isa.run(&cpu, 5));
    try testing.expectEqual(@as(u8, 5), cpu.a);
}

test "step and info table" {
    var cpu: Toy = .{};
    const r = try Toy.isa.step(&cpu, 0x01);
    try testing.expectEqual(@as(u8, 1), r.cycles);
    try testing.expectEqualStrings("add", Toy.isa.info[0x03].?.name);
    try testing.expectEqual(Toy.Mode.zp, Toy.isa.info[0x03].?.mode);
    try testing.expect(Toy.isa.info[0x00] == null);
}
