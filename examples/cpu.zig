//! A tiny accumulator CPU (`emu.dispatch.Isa`) wired through `emu.bus.Bus` to
//! RAM plus a memory-mapped output port. Runs a loop summing 1..10.
const std = @import("std");
const ts = @import("tsunami");
const dispatch = ts.emu.dispatch;
const bus = ts.emu.bus;

const Mode = enum { imp, imm, zp };

const Ram = struct {
    mem: [256]u8 = @splat(0),
    pub fn read(self: *Ram, addr: u16) u8 {
        return self.mem[addr & 0xFF];
    }
    pub fn write(self: *Ram, addr: u16, value: u8) void {
        self.mem[addr & 0xFF] = value;
    }
};

/// Records every byte written to it, in order — a fake "output port".
const OutPort = struct {
    log: [16]u8 = undefined,
    count: usize = 0,
    pub fn read(self: *OutPort, _: u16) u8 {
        return if (self.count == 0) 0 else self.log[self.count - 1];
    }
    pub fn write(self: *OutPort, _: u16, value: u8) void {
        self.log[self.count] = value;
        self.count += 1;
    }
};

const Cpu = struct {
    bus: *bus.Bus(2),
    pc: u16 = 0,
    a: u8 = 0,

    pub fn fetch(cpu: *Cpu) u8 {
        defer cpu.pc += 1;
        return cpu.bus.read8(cpu.pc);
    }
    fn operand(cpu: *Cpu, comptime mode: Mode) u8 {
        return switch (mode) {
            .imp => unreachable,
            .imm => cpu.fetch(),
            .zp => cpu.bus.read8(cpu.fetch()),
        };
    }

    pub fn lda(cpu: *Cpu, comptime mode: Mode) void {
        cpu.a = cpu.operand(mode);
    }
    pub fn add(cpu: *Cpu, comptime mode: Mode) void {
        cpu.a +%= cpu.operand(mode);
    }
    pub fn sta(cpu: *Cpu, comptime mode: Mode) void {
        std.debug.assert(mode == .zp);
        cpu.bus.write8(cpu.fetch(), cpu.a);
    }
    /// Decrement the zero-page counter at the fetched address; if the result
    /// is nonzero, add the signed (two's-complement) byte that follows to pc.
    pub fn dnz(cpu: *Cpu, comptime mode: Mode) void {
        std.debug.assert(mode == .zp);
        const addr = cpu.fetch();
        const v = cpu.bus.read8(addr) -% 1;
        cpu.bus.write8(addr, v);
        const rel: i8 = @bitCast(cpu.fetch());
        if (v != 0) cpu.pc = @bitCast(@as(i16, @bitCast(cpu.pc)) + rel);
    }
    pub fn out(cpu: *Cpu, comptime mode: Mode) void {
        std.debug.assert(mode == .zp);
        cpu.bus.write8(cpu.fetch(), cpu.a);
    }
    pub fn hlt(_: *Cpu, comptime _: Mode) dispatch.Stop {
        return .halt;
    }

    const isa = dispatch.Isa(Cpu, Mode, .{
        .{ .code = 0x01, .op = .lda, .mode = .imm, .cycles = 2 },
        .{ .code = 0x06, .op = .lda, .mode = .zp, .cycles = 3 },
        .{ .code = 0x02, .op = .add, .mode = .zp, .cycles = 3 },
        .{ .code = 0x03, .op = .sta, .mode = .zp, .cycles = 3 },
        .{ .code = 0x04, .op = .dnz, .mode = .zp, .cycles = 4 },
        .{ .code = 0x05, .op = .out, .mode = .zp, .cycles = 2 },
        .{ .code = 0xFF, .op = .hlt, .mode = .imp, .cycles = 1 },
    });
};

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    var ram: Ram = .{};
    var port: OutPort = .{};
    var b: bus.Bus(2) = .init();
    b.map(0x00, 0xEF, bus.Device.of(Ram, &ram));
    b.map(0xF0, 0xFF, bus.Device.of(OutPort, &port));

    // Zero page layout: [0x80]=counter (10), [0x81]=running sum.
    ram.mem[0x80] = 10;
    ram.mem[0x81] = 0;
    // sum = 0; while (counter != 0) { sum += counter; counter -= 1 }; out(sum); hlt
    const prog = [_]u8{
        0x01, 0, 0x03, 0x81, // LDA #0; STA sum
        0x06, 0x80, // top: LDA counter (zp)     (pc = 4)
        0x02, 0x81, // ADD sum (zp)
        0x03, 0x81, // STA sum
        0x04, 0x80, 247, // DNZ counter, rel=-9 -> back to pc=4 ("top")
        0x06, 0x81, // LDA sum (zp)
        0x05, 0xF0, // OUT sum -> port
        0xFF, // HLT
    };
    @memcpy(ram.mem[0..prog.len], &prog);

    var cpu: Cpu = .{ .bus = &b };
    const cycles = try Cpu.isa.run(&cpu, 10_000);

    try w.print("program: {d} bytes, {d} cycles used\n", .{ prog.len, cycles });
    try w.writeAll("disassembly:\n");
    var pc: usize = 0;
    while (pc < prog.len) {
        const op = prog[pc];
        const info = Cpu.isa.info[op] orelse {
            pc += 1;
            continue;
        };
        try w.print("  {x:0>2}: {s} ({s})\n", .{ pc, info.name, @tagName(info.mode) });
        pc += 1 + @as(usize, switch (info.mode) {
            .imp => 0,
            .imm, .zp => 1,
        });
        if (std.mem.eql(u8, info.name, "dnz")) pc += 1; // relative-jump byte
    }

    const sum: u32 = 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10;
    try w.print("output port received: {d} (want {d})\n", .{ port.log[port.count - 1], sum });
    if (port.count == 0 or port.log[port.count - 1] != sum) return error.ExampleFailed;
    try w.flush();
}
