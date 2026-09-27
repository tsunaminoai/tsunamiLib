const std = @import("std");
const contract = @import("../contract.zig");

pub const Device = struct {
    ptr: *anyopaque,
    vtable: struct {
        read: *const fn (*anyopaque, u16) u8,
        write: *const fn (*anyopaque, u16, u8) void,
    },

    /// Create a Device vtable from a type T with read(addr: u16) -> u8 and write(addr: u16, value: u8) -> void methods.
    pub fn of(comptime T: type, ptr: *T) Device {
        const ReadFn = *const fn (*T, u16) u8;
        const WriteFn = *const fn (*T, u16, u8) void;

        return .{
            .ptr = ptr,
            .vtable = .{
                .read = @ptrCast(@as(ReadFn, T.read)),
                .write = @ptrCast(@as(WriteFn, T.write)),
            },
        };
    }

    pub fn read(self: Device, addr: u16) u8 {
        return self.vtable.read(self.ptr, addr);
    }

    pub fn write(self: Device, addr: u16, value: u8) void {
        self.vtable.write(self.ptr, addr, value);
    }
};

pub const Region = struct {
    start: u16,
    end: u16,
    device: Device,
};

/// A memory bus with up to max_regions address ranges mapped to devices.
pub fn Bus(comptime max_regions: usize) type {
    return struct {
        regions: [max_regions]Region = undefined,
        region_count: usize = 0,

        pub fn init() @This() {
            return .{};
        }

        pub fn map(self: *@This(), start: u16, end: u16, device: Device) void {
            contract.require(self.region_count < max_regions, "bus: too many regions");
            contract.require(start <= end, "bus: invalid range");
            self.regions[self.region_count] = .{
                .start = start,
                .end = end,
                .device = device,
            };
            self.region_count += 1;
        }

        fn findRegion(self: @This(), addr: u16) ?*const Region {
            for (self.regions[0..self.region_count]) |*r| {
                if (addr >= r.start and addr <= r.end) {
                    return r;
                }
            }
            return null;
        }

        pub fn read8(self: @This(), addr: u16) u8 {
            if (self.findRegion(addr)) |region| {
                return region.device.read(addr);
            }
            return 0;
        }

        pub fn write8(self: *@This(), addr: u16, value: u8) void {
            if (self.findRegion(addr)) |region| {
                region.device.write(addr, value);
            }
        }

        pub fn read16le(self: @This(), addr: u16) u16 {
            const lo = self.read8(addr);
            const hi = self.read8(addr +| 1);
            return lo | (@as(u16, hi) << 8);
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const TestRam = struct {
    mem: [256]u8 = [_]u8{0} ** 256,

    pub fn read(self: *TestRam, addr: u16) u8 {
        return self.mem[addr & 0xFF];
    }

    pub fn write(self: *TestRam, addr: u16, value: u8) void {
        self.mem[addr & 0xFF] = value;
    }
};

const TestMirror = struct {
    base_addr: u16,
    actual_device: Device,

    pub fn read(self: *TestMirror, addr: u16) u8 {
        const actual_addr = self.base_addr + (addr & 0x7F);
        return self.actual_device.read(actual_addr);
    }

    pub fn write(self: *TestMirror, addr: u16, value: u8) void {
        const actual_addr = self.base_addr + (addr & 0x7F);
        self.actual_device.write(actual_addr, value);
    }
};

test Bus {
    var ram: TestRam = .{};
    var bus: Bus(2) = Bus(2).init();

    const device = Device.of(TestRam, &ram);
    bus.map(0x0000, 0x00FF, device);

    bus.write8(0x10, 42);
    const val = bus.read8(0x10);
    try std.testing.expectEqual(@as(u8, 42), val);
}

test "Bus read16le little-endian" {
    var ram: TestRam = .{};
    var bus: Bus(2) = Bus(2).init();

    const device = Device.of(TestRam, &ram);
    bus.map(0x0000, 0x00FF, device);

    bus.write8(0x20, 0x34);
    bus.write8(0x21, 0x12);
    const val = bus.read16le(0x20);
    try std.testing.expectEqual(@as(u16, 0x1234), val);
}

test "Bus returns 0 for unmapped address" {
    var bus: Bus(2) = Bus(2).init();
    const val = bus.read8(0xFF);
    try std.testing.expectEqual(@as(u8, 0), val);
}
