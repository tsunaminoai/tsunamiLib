const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;
const code_point = @import("code_point");
const ascii = std.ascii;

pub const Point = struct {
    x: usize = 0,
    y: usize = 0,
};

pub const Version = enum(u8) {
    invalid = 0,
    @"1" = 1,
    @"2" = 2,
    @"32" = 32,
    @"40" = 40,
    pub fn get_config(self: Version) Config {
        return switch (@intFromEnum(self)) {
            1...9 => .{ .bits = 148, .codewords = 20 },
            10...26 => .{ .bits = 156, .codewords = 20 },
            27...40 => .{ .bits = 156, .codewords = 20 },
            else => unreachable,
        };
    }
    pub fn asInt(self: Version) u8 {
        return @intFromEnum(self);
    }
    pub const Config = struct {
        ecc: enum { L, M, Q, H } = .L,
        bits: usize = 0,
        codewords: usize = 0,
    };
    pub fn getNumRawDataModules(self: Version) usize {
        const i: usize = @intCast(self.asInt());
        var res: usize = (16 * i + 128) * i + 64;
        if (i >= 2) {
            const numAlign = @divFloor(i, 7) + 2;
            res -= (25 * numAlign - 10) * numAlign - 55;
            if (i >= 7) res -= 36;
        }
        return res;
    }
    pub fn getNumDataCodewords(self: Version) usize {
        const cfg = self.get_config();
        return @intCast(@as(isize, @intCast(@divFloor(self.getNumRawDataModules(), 8))) -
            ECCCodeWordsPerBlock[@intFromEnum(cfg.ecc)][@intCast(self.asInt())] *
                NumECCBlocks[@intFromEnum(cfg.ecc)][@intCast(self.asInt())]);
    }
};

pub const Segment = struct {
    mode: Mode, // The segment mode is always a 4-bit field.
    count: u8, // The character count’s field width depends on the mode and version.
    numChars: usize, // nuumber of characters encoded in the bitfield
    data: []const u1,
    terminator: u4 = 0, // The terminator is normally four “0” bits, but fewer if the data codeword capacity is reached.
    version: Version,
    //The bit padding is between zero to seven “0” bits, to fill all unused bits in the last byte.
    // The byte padding consists of alternating (hexadecimal) EC and 11 until the capacity is reached.

    pub const Mode = enum(u4) {
        number = 0x1,
        alpha = 0x2,
        byte = 0x4,
        kanji = 0x8,
        eci = 0x7,
        const Numeric = [_]usize{ 10, 12, 14 };
        const Alpha = [_]usize{ 9, 11, 13 };
        const Byte = [_]usize{ 8, 16, 16 };
        const Kanji = [_]usize{ 8, 10, 12 };
        const Eci = [_]usize{ 0, 0, 0 };
        pub fn numCharCountBits(self: Mode, ver: Version) usize {
            return (switch (self) {})[@intCast(@divFloor(ver.asInt() + 7, 17))];
        }
        const UniTest = packed struct(u4) {
            is_num: bool = false,
            is_alpha: bool = false,
            is_byte: bool = false,
            is_kanji: bool = false,
            pub fn asInt(self: UniTest) u4 {
                return @bitCast(self);
            }
            pub fn orAll(self: *UniTest, other: UniTest) void {
                self.* = @bitCast(self.asInt() | other.asInt());
            }
        };
        pub fn analyze_unicode(input: []const u8) !Mode {
            var tests = UniTest{};
            var iter: code_point.Iterator = .init(input);
            while (iter.next()) |cp| {
                // if (!uc.utf8ValidCodepoint(cp)) return error.InvalidUnicode;
                var t = UniTest{
                    .is_num = isNumber(cp),
                    .is_alpha = isAlpha(cp),
                    .is_kanji = isKanji(cp),
                };
                if (t.asInt() == 0) t.is_byte = true;
                tests.orAll(t);
            }
            return if (tests.is_kanji) .kanji else if (tests.is_byte) .byte else if (tests.is_num) .number else .alpha;
        }
        fn isNumber(cp: code_point.CodePoint) bool {
            const zero = code_point.decode("0", 0).?;
            const nine = code_point.decode("0", 0).?;
            return zero.code <= cp.code and cp.code <= nine.code;
        }
        fn isAlpha(cp: code_point.CodePoint) bool {
            const char_set = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";

            return cp.code < 128 and std.mem.indexOf(u8, char_set, &.{@intCast(cp.code)}) != null;
        }

        fn isKanji(cp: code_point.CodePoint) bool {
            if (cp.code >= 0x10000) return false;
            const trunc = @as(u16, @truncate(cp.code >> 2));
            const p = std.fmt.parseInt(
                u16,
                KanjiSet[trunc .. trunc + 1],
                16,
            ) catch return false;
            return 0 != (p >> @intCast(trunc & @as(u3, 3))) & @as(u3, 1);
        }
    };
    pub fn init(m: Mode, ver: Version, bit_data: []const u1, num_chars: usize) !Segment {
        if (num_chars == 0) return error.InvalidNumberofChars;
        const s = Segment{
            .mode = m,
            .numChars = num_chars,
            .data = bit_data,
            .version = ver,
            .count = @intCast(ver.get_config().bits),
        };
        return s;
    }
    fn getTotalBits(segments: []const Segment, ver: Version) ?usize {
        var res: usize = 0;
        for (segments) |seg| {
            const ccbits = seg.mode.numCharCountBits(ver);
            if (seg.numChars >= (@as(usize, 1) << @intCast(ccbits)))
                return null;
            res += 4 + ccbits + (if (seg.data) |d| d.len else 0);
        }
        return res;
    }
};
test "Segments" {
    const input = Tests.hello_world;
    _ = input; // autofix
    try tst.expectEqual(.byte, try Segment.Mode.analyze_unicode(Tests.hello_world));
    try tst.expectEqual(.kanji, Segment.Mode.analyze_unicode(Tests.kanji));
    try tst.expectEqual(.number, try Segment.Mode.analyze_unicode(Tests.numeric));
    try tst.expectEqual(.byte, try Segment.Mode.analyze_unicode(Tests.utf8));
    try tst.expectEqual(.alpha, try Segment.Mode.analyze_unicode(Tests.alpha));
    const s = try Segment.init(.byte, .@"2", &[_]u1{ 1, 0, 1, 0, 1 }, 1);
    _ = s; // autofix
    // try tst.expectEqual(1, s.count);
}

pub const Block = struct {};

pub const Codeword = struct {
    value: u8,

    preInterleaveIndex: isize = -1,
    blockIndex: isize = -1,
    indexInBlock: isize = -1,
    postInterleaveIndex: isize = -1,

    pub fn init(v: u8) !Codeword {
        if (v < 0 or v > 255) return error.InvalidValue;
        return .{
            .value = v,
        };
    }
};
const DataCodeWord = struct {
    preEccIndex: isize = -1,
    cw: Codeword,

    pub fn init(value: u8) !DataCodeWord {
        return .{ .cw = try .init(value) };
    }
};

pub const Module = union(enum) {
    unfilled: Base,
    function: Function,
    code_word: CodeWord,
    remainder: Base,
    mask: Mask,
    filled: Filled,

    pub const Base = struct {
        color: bool = false,
    };
    pub const Unfilled = Base{};
    pub const CodeWord = Base;
    pub const Remainder = Base{ .color = false };
    pub const Mask = Base;

    pub const Filled = struct {
        is_new: bool,
        color: bool,
    };

    pub const Function = struct {
        kind: Kind = .separator,
        color: bool,
        pub const Kind = enum {
            finder,
            separator,
            timing,
            alignment,
            format_info,
            version_info,
            dark,
        };
    };
};

inline fn getBit(x: anytype, i: usize) u1 {
    return @intCast((x >> @intCast(i)) & 1);
}
test "getbit" {
    try tst.expectEqual(1, getBit(0b00100, 2));
}

pub const QRCode = struct {
    side_len: usize,
    modules: Array(Array(Module)), // side x side
    version: Version,
    ecc_level: ECCLevel,

    allocator: Allocator,

    pub fn init(a: Allocator, size: usize, ecc: ECCLevel, version: Version) !QRCode {
        var self = QRCode{
            .allocator = a,
            .side_len = size,
            .ecc_level = ecc,
            .version = version,
            .modules = Array(Array(Module)).init(a),
        };
        for (0..size) |y| {
            try self.modules.append(Array(Module).init(a));
            try self.modules.items[y].appendNTimes(.{ .unfilled = .{} }, size);
        }

        return self;
    }

    pub fn deinit(self: *QRCode) void {
        for (self.modules.items) |row| {
            row.deinit();
        }
        self.modules.deinit();
    }

    pub fn clearNewFlags(self: *QRCode) void {
        for (self.modules.items) |*column| {
            for (column.items) |*mod| {
                switch (mod.*) {
                    .filled => |*f| f.is_new = false,
                    else => continue,
                }
            }
        }
    }

    pub fn makeZigZagScan(self: QRCode) ![]Point {
        var res = Array(Point).init(self.allocator);
        defer res.deinit();
        var right: isize = self.side_len;
        while (right >= 1) : (right -= 2) {
            if (right == 6) right = 5;
            var vert: usize = 0;
            while (vert < self.side_len) : (vert += 1) {
                for (0..2) |j| {
                    const x = right - j;
                    const up = ((right + 1) & 2) == 0;
                    const y = if (up) self.side_len - 1 - vert else vert;
                    if (self.modules.items[x].items[y] == .function)
                        try res.append(.{ x, y });
                }
            }
        }
        return try res.toOwnedSlice();
    }
    pub fn applyMask(self: *QRCode, mask: QRCode) !void {
        for (0..self.side_len) |x| {
            for (0..self.side_len) |y| {
                const a = mask.modules.items[x].items[y];
                const b = self.modules.items[x].items[y];
                if (a == .mask and b == .filled)
                    b.color = b.color != a.color;
            }
        }
    }
    pub fn drawCodewords(self: QRCode, codewords: []const Codeword, zig_zag: []const Point) !void {
        if (codewords.len != @divFloor(self.version.getNumRawDataModules(), 8)) {
            return error.InvalidArgument;
        }
        for (zig_zag, 0..) |z, i| {
            if (i < codewords.len * 8) {
                const cw = codewords[i >> 3];
                self.modules.items[z.x].items[z.y] = .{
                    .code_word = .init(getBit(cw.value, 7 - (i & 7))),
                };
            } else self.modules.items[z.x].items[z.y] = .{ .remainder = .{} };
        }
    }
    pub fn interleaveBlocks(self: QRCode, blocks: []const []const Codeword) ![]ECCCodeWord {
        const block_len = ECCCodeWordsPerBlock[@intFromEnum(self.ecc_level)][@intFromEnum(self.version)];
        _ = block_len; // autofix
        const short_len = blocks[0].len;
        _ = short_len; // autofix

    }

    pub fn splitIntoBlocks(self: QRCode, data: []const DataCodeWord) !Array([]const DataCodeWord) {
        const numBlocks = NumECCBlocks[@intFromEnum(self.ecc_level)][@intFromEnum(self.version)];
        const eccBlockLen = ECCCodeWordsPerBlock[@intFromEnum(self.ecc_level)][@intFromEnum(self.version)];
        const rawCodeWords = @divFloor(self.version.getNumRawDataModules(), 8);
        const numShortBlocks = numBlocks - @mod(rawCodeWords, numBlocks);
        const shortBlockLen = @divFloor(rawCodeWords, numBlocks);

        var res = Array([]const DataCodeWord).init(self.allocator);
        errdefer res.deinit();

        var off: usize = 0;
        for (0..numBlocks) |bi| {
            const end = off + shortBlockLen - eccBlockLen + if (bi < numShortBlocks) 0 else 1;
            const block = data[off..end];
            for (block, 0..) |*dcw, i| {
                dcw.cw.blockIndex = bi;
                dcw.cw.indexInBlock = i;
            }
            try res.append(block);
            off = end;
        }
        return res;
    }

    pub fn drawVersionInfo(self: *QRCode) !void {
        if (@intFromEnum(self.version) < 7) return;

        var rem: usize = @intFromEnum(self.version);
        for (0..12) |_| {
            rem = (rem << 1) ^ ((rem) >> 11) * 0x1F25;
        }
        const bits: usize = @as(usize, @intCast(self.version.asInt())) << 12 | rem;
        if (bits >> 18 != 0) return error.WrongBits;

        for (0..18) |i| {
            const color = getBit(bits, i) == 1;
            const a = self.side_len - 11 + @mod(i, 3);
            const b = @divFloor(i, 3);
            self.modules.items[a].items[b] = .{ .function = .{ .kind = .version_info, .color = color } };
            self.modules.items[b].items[a] = .{ .function = .{ .kind = .version_info, .color = color } };
        }
    }

    pub fn drawTimingPatterns(self: *QRCode) void {
        for (0..self.side_len) |i| {
            self.modules.items[6].items[i] = .{ .function = .{ .kind = .timing, .color = @mod(i, 2) == 0 } };
            self.modules.items[i].items[6] = .{ .function = .{ .kind = .timing, .color = @mod(i, 2) == 0 } };
        }
    }
    pub fn drawFinderPatterns(self: *QRCode) void {
        const centers = [_][]const usize{
            &.{ 3, 3 },
            &.{ self.side_len - 4, 3 },
            &.{ 3, self.side_len - 4 },
        };
        for (centers) |c| {
            var dy: isize = -4;
            while (dy <= 4) : (dy += 1) {
                var dx: isize = -4;
                while (dx <= 4) : (dx += 1) {
                    const dist = @max(@abs(dx), @abs(dy));
                    const x = @as(isize, @intCast(c[0])) + dx;
                    const y = @as(isize, @intCast(c[1])) + dy;
                    if (0 <= x and x < self.side_len and 0 <= y and y < self.side_len) {
                        self.modules.items[@intCast(x)].items[@intCast(y)] = .{
                            .function = .{
                                .kind = if (dist <= 3) .finder else .separator,
                                .color = dist != 2 and dist != 4,
                            },
                        };
                    }
                }
            }
        }
    }

    pub fn drawFormatbits(self: *QRCode, mask: anytype) !void {
        var bits: usize = 0;
        if (mask != -1) {
            const data = self.ecc_level.getFormatBits() << 3 | mask;
            var rem = data;
            for (0..10) |_|
                rem = (rem << 1) ^ ((rem >> 9) * 0x0537);

            bits = (data << 10 | rem) ^ 0x5412;
        }
        if (bits >> 15 != 0) return error.Assertion;
        for (0..5) |i|
            self.setFmtInfoModule(8, i, getBit(bits, i) == 1);
    }
    fn setFmtInfoModule(self: *QRCode, x: usize, y: usize, color: bool) void {
        self.modules.items[x].items[y] = .{ .function = .{ .kind = .format_info, .color = color } };
    }
    pub fn drawAlignmentPatterns(self: *QRCode) !void {
        if (self.version == .@"1") return;
        var positions = Array(isize).init(self.allocator);
        defer positions.deinit();
        try positions.append(6);
        const numAlign = @divFloor(self.version.asInt(), 7) + 2;
        const step: isize = if (self.version == .@"32")
            26
        else
            @intCast(try math.divCeil(
                usize,
                self.side_len - 13,
                numAlign * 2 - 2,
            ) * 2);
        var pos: isize = @as(isize, @intCast(self.side_len)) - 7;
        while (positions.items.len < numAlign) : (pos -= step) {
            //todo: splice?
            // positions.splice(1, 0, pos);
            try positions.appendSlice(&.{ 1, 0, pos });
        }
        for (positions.items, 0..) |cx, i| {
            for (positions.items, 0..) |cy, j| {
                if ((i == 0 and j == 0) or (i == 0 and j == numAlign - 1) or (i == numAlign - 1 and j == 0))
                    return;

                for (0..4) |d_x| {
                    const dx: isize = 2 - @as(isize, @intCast(d_x));
                    for (0..4) |d_y| {
                        const dy: isize = 2 - @as(isize, @intCast(d_y));
                        self.modules.items[@intCast(cx + dx)].items[@intCast(cy + dy)] = .{ .function = .{ .kind = .alignment, .color = @max(@abs(dx), @abs(dy)) != 1 } };
                    }
                }
            }
        }
    }
    pub fn makeMask(self: QRCode, mask: usize) !QRCode {
        var res = try QRCode.init(
            self.allocator,
            self.side_len,
            self.ecc_level,
            self.version,
        );
        errdefer res.deinit();
        for (0..self.side_len) |x| {
            for (0..self.side_len) |y| {
                const invert = switch (mask) {
                    0 => @mod(x + y, 2) == 0,
                    1 => @mod(y, 2) == 0,
                    2 => @mod(x, 3) == 0,
                    3 => @mod(x + y, 3),
                    4 => @divFloor(x, 3) + @mod(@divFloor(y, 2), 2) == 0,
                    5 => x * @mod(y, 2) + x * @mod(y, 3) == 0,
                    6 => @mod(x * @mod(y, 2) + x * @mod(y, 3), 2) == 0,
                    7 => @mod(@mod(x + y, 2) + x * @mod(y, 3), 2) == 0,
                    else => return error.InvalidMask,
                };
                if (self.modules.items[x].items[y] == .function) {
                    self.modules.items[x].items[y] = .{ .mask = .{ .color = invert } };
                }
            }
        }
        return res;
    }
    pub fn computePenalties(self: QRCode) !Penalty.Info {
        const penalties: [4]usize = [_]usize{0} ** 4;
        var info = Penalty.Info.init(self.allocator);
        errdefer info.deinit();
        const colors = Array([]bool).init(self.allocator);
        defer colors.deinit();

        for (self.modules.items) |col| {
            for (col.items) |cell| {
                try colors.append(cell == .filled and cell.color);
            }
        }
        for (0..self.side_len) |y| {
            var run_color = false;
            var runX: usize = 0;
            var runHist = FinderPenalty.init(
                self.allocator,
                self.side_len,
                .horiz,
                y,
                info.horizFalseFinders,
            );
            errdefer runHist.deinit();
            for (0..self.side_len) |x| {
                if (colors.items[x][y] == run_color) {
                    runX += 1;
                    if (runX == 5) {
                        penalties[0] += Penalty.n1.asInt(usize);
                        try info.horizRuns.append(.init(x + 1 - runX, y, runX));
                    } else if (runX > 5) {
                        penalties[0] += 1;
                        info.horizRuns.items[info.horizRuns.items.len - 1] = LinearRun.init(
                            x + 1 - runX,
                            y,
                            runX,
                        );
                    }
                } else {
                    try runHist.addHistory(runX);
                    if (!run_color)
                        penalties[2] ++ (try runHist.countAndAddPatterns()) * Penalty.n3.asInt(usize);

                    run_color = colors.items[x][y];
                    runX += 1;
                }
            }
        }

        for (0..self.side_len) |x| {
            var run_color = false;
            var runY: usize = 0;
            var runHist = FinderPenalty.init(
                self.allocator,
                self.side_len,
                .vert,
                x,
                info.horizFalseFinders,
            );
            errdefer runHist.deinit();
            for (0..self.side_len) |y| {
                if (colors.items[x][y] == run_color) {
                    runY += 1;
                    if (runY == 5) {
                        penalties[0] += Penalty.n1.asInt(usize);
                        try info.vertRuns.append(.init(x, y + 1 - runY, runY));
                    } else if (runY > 5) {
                        penalties[0] += 1;
                        info.vertRuns.items[info.vertRuns.items.len - 1] = LinearRun.init(
                            x,
                            y + 1 - runY,
                            runY,
                        );
                    }
                } else {
                    try runHist.addHistory(runY);
                    if (!run_color)
                        penalties[2] ++ (try runHist.countAndAddPatterns()) * Penalty.n3.asInt(usize);

                    run_color = colors.items[x][y];
                    runY += 1;
                }
            }
        }

        for (0..self.side_len - 1) |x| {
            for (0..self.side_len - 1) |y| {
                const c = colors.items[x][y];
                if (c == colors.items[x + 1][y] and c == colors[x][y + 1] and c == colors[x + 1][y + 1]) {
                    penalties[1] += Penalty.n2.asInt(usize);
                    try info.twoByTwoBoxes.append(.{ .x = x, .y = y });
                }
            }
        }

        var dark: usize = 0;
        for (colors.items) |col| {
            dark = @reduce(.Add, col);
        }
        const total = self.side_len * self.side_len;
        var k: usize = 0;
        while (@abs(dark * 20 - total * 10) > (k + 1) * total)
            k += 1;
        penalties[3] += k * Penalty.n4.asInt(usize);
        @memcpy(&info.penalties, &penalties);
        return info;
    }
    pub fn format(self: QRCode, comptime fmt: []const u8, options: anytype, writer: anytype) !void {
        _ = fmt; // autofix
        _ = options; // autofix
        const space = "  ";
        const block = "██";
        for (self.modules.items) |row| {
            for (row.items) |cell| {
                const out = blk: switch (cell) {
                    .filled => |c| break :blk if (c.color) block else space,
                    .function => |f| {
                        break :blk switch (f.kind) {
                            .separator => space,
                            inline else => if (f.color) block else space,
                        };
                    },
                    .mask => |c| break :blk if (c.color) block else space,
                    inline else => break :blk space,
                };
                try writer.writeAll(out);
            }
            try writer.writeAll("\n");
        }
    }
};

// from https://www.nayuki.io/page/creating-a-qr-code-step-by-step
const Tests = .{
    .hello_world = "Hello, world! 123",
    .kanji = "「魔法少女まどか☆マギカ」って、　ИАИ　ｄｅｓｕ　κα？",
    .numeric = "31415926535897932384626433832795028841971693993",
    .utf8 = "aЉ윇😱",
    .alpha = "PROJECT NAYUKI",
};

test "version" {
    const v2 = Version.@"2";
    try tst.expectEqual(34, v2.getNumDataCodewords());
}

test "hello world" {
    const input = Tests.hello_world;
    _ = input; // autofix
    // std.debug.print("{}\n", .{cfg});
    const expect_codewords = [_]u8{ 0x41, 0x14, 0x86, 0x56, 0xC6, 0xC6, 0xF2, 0xC2, 0x07, 0x76, 0xF7, 0x26, 0xC6, 0x42, 0x12, 0x03, 0x13, 0x23, 0x30, 0x85, 0xA9, 0x5E, 0x07, 0x0A, 0x36, 0xC9 };
    _ = expect_codewords; // autofix

    var q = try QRCode.init(tst.allocator, 25, .low, .@"2");
    defer q.deinit();
    q.clearNewFlags();
    q.drawTimingPatterns();
    q.drawFinderPatterns();
    try q.drawAlignmentPatterns();
    try q.drawVersionInfo();
    // try q.drawFormatbits(2);

    std.debug.print("{}\n", .{q});
}

const KanjiSet = "0000000000000000000000010000000000000000C811350000000800000008000000000000000000000000000000000000000000000000000000000000000000" ++
    "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000EFFFBF30EFFFBF30000000000000" ++
    "2000FFFFFFFFFFFFFFFF200000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" ++
    "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" ** 13 ++
    "000016333600D080000000000000000000000000000000000000000000000000800000000080000000000000000000000000F000000000000000410000000000" ++
    "D890404618F10302000040003CC00000CC0000000200000000000000000000000000400000000000000000000000000000000000000000000000000000000000" ++
    "0000000000000000000000000000000000000000000000000000000000000000F0099993939999994080000000000000000000003000C0030C8C000000080000" ++
    "060000000000000050000000004A0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" ++
    "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" ** 4 ++
    "FEFFF30100000000EFFFFFFFFFFFFFFFFFFFF087EFFFFFFFFFFFFFFFFFFFF7870000000000000000000000000000000000000000000000000000000000000000" ++
    "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" ** 14 ++
    "B8F63F34244264B9C28E0E3E4000A00456F563BD779794407DCE0F3E83065C80206E3043000815538C0EBAE700289689849A2492308E0608C14439DAA30C8654" ++
    "06AA6568A7F304208838164102014712120220700003CB04426A26448A0602A071204758204048C9BFB7514142F72E11566BFE2057F1FF0207A304833C623676" ++
    "9DD429020B649CF089CB05848368F30A8832618890E32325AA224A3EDD00C27C661A1E62B048A0F8BE72E955142CBB984100045816369480C0F70DA8E3FFFC50" ++
    "A1FF308A14A704B7547420080050BE83158D50004399C01779300010663640420D081500000CA03417098C038000008599E0007F08F514000B00014981000826" ++
    "04200D90002865104005108001D101501C48010052040501F014A8D49004D06A91BAC4190C121890584C30002560000840B08000D14090009484C50990000961" ++
    "56C002222148334230C0697440A0522482008809009480F42A41AA3D038D78E3406816F14AE76814093C3B505A758112E14284A2821140A404A0B16106D00488" ++
    "A0202059122806013420004044410008000040C00000000760A11C00A42000C000A104004041540492003BDB87A0B2509ABB0AFBC7049738CF21D18E6FB4965C" ++
    "6FFEA440511220FF36DEB204330D24200001310020B1AC950A000020307A14C208842FF840200000008550010029840049811002400508430023C486AE94EB86" ++
    "C48124E2028A9C129B050B08E100C7FFA9480E411C820E10E07894CAF031BDDDA1EBBF980E2A2A152055AC2364E3B829FBD1F8343076812382030C331180329A" ++
    "000C56A33EF82040E4C25268D3FB1A00A1A34DC89C60C7130E00A059B810BDE0B43E02C82811010F49D7827ACA9CBF044844356009A544448CF3100084004D5F" ++
    "107775CE244CD19838B682949014242DD160EF95008122A34E7BF9B3300FAE0C683120280898004E002B1A0108B44CC0903D4498FAF14384952854C0A0240540" ++
    "040A8C010413054440040010082804508010C24403A650A16A024150FC0965461200001381C90FBC021A2E36C4015B10C83538A92B8B1823A78948A07E3320C0" ++
    "CC4D81091A1A0709E1A8400E4D3C1540A9342C1244840135292004631420DB3F90BA0F8E0CD72D5A242CB42DF34AFA0D0AA11A4374288D30254CB156492DA38C" ++
    "C1008C0460E04133F416B12B88000D0CA20A898A5C1AB66105E24B58B80C4060339F40E1E650152A0040836770CE8B37604423811804618CA8C79036089240AA" ++
    "42C1C9ACE0E40672099000386400401811061801D0458090E000A0CC005000C00340440AB8004876882591A8E56881B895E2061401C8EBC91686C19800898000" ++
    "0018A9808100091470041A4E5050D046E013D4E06084A0FF23618AA2E258B000008148AC02E0C6962300006185650930021582A1000842111E81623425D5AAE0" ++
    "0AF082EAB7AF005480460498088C440C5009141B42484C4243A1A3060009491C6428A300D081601C22000199050E115175042800A140A020F4000398318DA444" ++
    "20A822DE0C015004000120108088101300644020000F80700098002A000020220020016124000401002506204F250002015803280011202480345B081E0702A9" ++
    "04021080005356CF1C9140BA682041267800440058094420C50458A07023083300400C8B02EC0D0C030C0800805052D009A004000020C0805056000412462014" ++
    "862000004200C748200002ED91689404808000044800100200480101DC247C108307A25D8691F8D105EB21E35FE29D184CEC21428280E237CA4243B4C020D14D" ++
    "20A20008790011804C114411687154D79D94946000041978C4524C8DAB44419429B1008C17200851180000C0A690002C00842004120394AB080208C1CA2E8001" ++
    "400143001E00414802000002008941012C07AA408868024526C03140081901022804602004C1004538309E4E52120848334E00020C44906E30A06218AD211080" ++
    "109609791004688FD42E1800E0A0156AA110CE18006C14488BDAC26BF64A147845D820B41611862006BB75020A0533400C8A4B7B204221103DA9000217228C00" ++
    "1802E908A8C0081E900B151813018204E0A25A986B96E0265244441D580845D457C21BF1708DD268C78D1484E414E622002880E9C08F73DE08C8625731394180" ++
    "23E0408CE4846AE6A4C207660C6210ABC03DD58100000000000000000000000000000000000004500207331818F45A30CE550146813C44322641430034A090A1" ++
    "B7815A312010848A0440445C6018DD2E0FA184D2626B6140850504E6230821134E7000C08028A0240484B30806843178D05882439130925E5432A0789A5281C5" ++
    "6A775C9D00E58E301800007A45DC6C140082402A068BA0B20200E9ADAE80A1E0C7240C11809867301268501000008A810A64844C50D022022245841910A87982" ++
    "898780005061221304240880E4AF2A6041801129200240B925E23460000501091400AB5824030202B0F40A5080720804195039A105FD0060202A1203036008E4" ++
    "4CC08C40400A10080006134D088020A000C1820081E814000DA601AC012F00B4D47260510220098800A58A289454051840020008408880C21D80500010084CA4" ++
    "020E2600E80024A05503C8A6E0905A0E60924C2441843B08E308033B2010C1374941D00602C00490A103040C154A490CACD88C502C69C04A100040438C000110" ++
    "D0559C9A8242A5C124107384D4A7F0014B23A254B7135029498B44C57D86A85479051DE234918024202D04D9048979029045D460000000000000000000000000" ++
    "00000000000008482455124060C100714267011678FFDD9BF956A0C5D54200C30BDA950000000000000000000D82B9002240002080108044408904CAA0D88209" ++
    "0078100E0040130049711228910045012BC2A12020C9300406D34088C08000978261C3AB046880BC47270809E10000000000008D881E78C94304214046EA1972" ++
    "B68EBF6EF80466356AEEF735B23E4E5BF9682000845822102044008061120A02400040200002500000E74510C261CA1048A2580141C803503CBF349BAC000609" ++
    "000623040021090803B018C4450020049200A6D100020820000840000162C05104081070D49D42C0018205811005020500010D400807102080103C1223100000" ++
    "88009170208006502100E0C450401A0F2000000000000000000000000000000000000000000000000000000000000800D8E8A530DB1240A58843071162000000" ++
    "00000001205C4088518B108FC741DE5206DE0BB198507DB13FA726A1C0D05CA01D5EA425094050364530442575B22161278A101194928100849080010006C688" ++
    "E619F85021030993048F03940888B10000000000005824008500008940AE41078261D1163115000642A17A000000000000000C30021781012710729A40066098" ++
    "220CC02000901804D2020AC843E0000000000000001210111108A11CC4CE298004000058CA7C6081E30E2150000801008004EC0810D6012014686580E1107200" ++
    "0573D380230E50E40C104840180004100000000000000000000000000AA195008C34428884D1008C25103027310940400828004001A841D065088020040A4072" ++
    "000000C400000000000000000000023A2091EA0A066200FD010F51B712180DA3081482003001008400CC4108FC414C0000020203100000000000000000000000" ++
    "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" ** 47 ++
    "0000000000000000000000000000000000000000000000000000000000000000A7FDFFFFFFFFFFFEFFFFFFF30000000000000000000000000000000082000000";

comptime {
    std.debug.assert(KanjiSet.len == 16384);
}

pub const ECCLevel = union(enum) {
    low,
    med,
    qrt,
    hi,
    pub fn getFormatBits(self: ECCLevel) usize {
        _ = self; // autofix
        return 34;
    }
};
pub const ECCCodeWord = DataCodeWord;
pub const ECCCodeWordsPerBlock = [_][]const i16{
    &.{ -1, 7, 10, 15, 20, 26, 18, 20, 24, 30, 18, 20, 24, 26, 30, 22, 24, 28, 30, 28, 28, 28, 28, 30, 30, 26, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 },
    &.{ -1, 10, 16, 26, 18, 24, 16, 18, 22, 22, 26, 30, 22, 22, 24, 24, 28, 28, 26, 26, 26, 26, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28 },
    &.{ -1, 13, 22, 18, 26, 18, 24, 18, 22, 20, 24, 28, 26, 24, 20, 30, 24, 28, 28, 26, 30, 28, 30, 30, 30, 30, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 },
    &.{ -1, 17, 28, 22, 16, 22, 28, 26, 26, 24, 28, 24, 28, 22, 24, 24, 30, 28, 28, 26, 28, 30, 24, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30 },
};
pub const NumECCBlocks = [_][]const i16{
    &.{ -1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 7, 8, 8, 9, 9, 10, 12, 12, 12, 13, 14, 15, 16, 17, 18, 19, 19, 20, 21, 22, 24, 25 },
    &.{ -1, 1, 1, 1, 2, 2, 4, 4, 4, 5, 5, 5, 8, 9, 9, 10, 10, 11, 13, 14, 16, 17, 17, 18, 20, 21, 23, 25, 26, 28, 29, 31, 33, 35, 37, 38, 40, 43, 45, 47, 49 },
    &.{ -1, 1, 1, 2, 2, 4, 4, 6, 6, 8, 8, 8, 10, 12, 16, 12, 17, 16, 18, 21, 20, 23, 23, 25, 27, 29, 34, 34, 35, 38, 40, 43, 45, 48, 51, 53, 56, 59, 62, 65, 68 },
    &.{ -1, 1, 1, 2, 4, 4, 4, 5, 6, 8, 8, 11, 11, 16, 16, 18, 16, 19, 21, 25, 25, 25, 34, 30, 32, 35, 37, 40, 42, 45, 48, 51, 54, 57, 60, 63, 66, 70, 74, 77, 81 },
};
pub const ReedSolomonGenerator = struct {
    coeffs: Array(u8),
    allocator: Allocator,
    pub fn init(a: Allocator, degree: usize) ReedSolomonGenerator {
        if (degree < 1 or degree > 255) return error.DegreeOutOfRange;
        var self = ReedSolomonGenerator{
            .coeffs = Array(u8).init(a),
            .allocator = a,
        };
        for (0..degree - 1) |_|
            try self.coeffs.append(0);

        try self.coeffs.append(1);

        var root: usize = 1;
        for (0..degree) |_| {
            for (0..self.coeffs.items.len) |j| {
                self.coeffs.items[j] = self.multiply(self.coeffs.items[j], root);
                if (j + 1 < self.coeffs.items.len)
                    self.coeffs.items[j] ^= self.coeffs.items[j + 1];
            }
            root = self.multiply(root, 0x2);
        }

        return self;
    }
    pub fn deinit(self: ReedSolomonGenerator) void {
        self.coeffs.deinit();
    }

    fn multiply(self: ReedSolomonGenerator, x: u8, y: u8) !u8 {
        _ = self; // autofix
        if (x >> 8 != 0 or y >> 8 != 0) return error.ByteOutOfRange;

        var z: u8 = 0;
        var i: isize = 7;
        while (i >= 0) : (i -= 1) {
            z = (z << 1) ^ ((z >> 7) * 0x11D);
            z ^= ((y >> i) & 1) * x;
        }
        if (z >> 8 != 0) return error.GotZero;
        return z;
    }
    pub fn getRemainder(self: *ReedSolomonGenerator, data_bytes: []const u8) ![]u8 {
        var res = try self.coeffs.clone();
        @memset(&res.items, 0);

        for (data_bytes) |byte| {
            const factor = byte ^ (res.orderedRemove(0));
            try res.append(0);
            for (self.coeffs.items, 0..) |c, i| {
                res.items[i] ^= self.multiply(c, factor);
            }
        }

        return try res.toOwnedSlice();
    }
};

pub const Penalty = enum(u8) {
    n1 = 3,
    n2 = 3,
    n3 = 40,
    n4 = 10,
    pub fn asInt(self: Penalty, comptime T: type) T {
        return @intCast(self);
    }
    pub const Points = struct { usize, usize, usize, usize };
    pub const Info = struct {
        horizRuns: Array(LinearRun),
        vertRuns: Array(LinearRun),
        twoByTwoBoxes: Array(Point),
        horizFalseFinders: Array(LinearRun),
        vertFalseFinders: Array(LinearRun),
        numDarkModules: usize = 0,
        penaltyPoints: Array(Points),
        pub fn init(a: Allocator) Info {
            return .{
                .horzRuns = Array(Point).init(a),
                .vertRuns = Array(Point).init(a),
                .twoByTwos = Array(Point).init(a),
                .horzFinders = Array(Point).init(a),
                .vertFinders = Array(Point).init(a),
                .dark = Array(Point).init(a),
            };
        }
        pub fn deinit(self: Info) void {
            self.horizRuns.deinit();
            self.vertRuns.deinit();
            self.twoByTwoBoxes.deinit();
            self.horizFalseFinders.deinit();
            self.vertFalseFinders.deinit();
            self.penaltyPoints.deinit();
        }
    };
};
pub const LinearRun = struct {
    startX: usize = 0,
    startY: usize = 0,
    runLen: usize = 0,
    pub fn init(x: usize, y: usize, len: usize) LinearRun {
        return .{
            .startX = x,
            .startY = y,
            .runLen = len,
        };
    }
};
pub const FinderPenalty = struct {
    runHistory: Array(usize),
    runEndPositions: Array(usize),
    position: usize = 0,
    padding: usize,
    direction: Dir,
    qr_size: usize,

    outer: usize,
    finders: Array(LinearRun),

    const Dir = enum { horiz, vert };

    pub fn init(
        a: Allocator,
        qr_size: usize,
        dir: Dir,
        outer_position: usize,
        finders: Array(LinearRun),
    ) FinderPenalty {
        return .{
            .runHistory = Array(usize).init(a),
            .runEndPositions = Array(usize).init(a),
            .outer = outer_position,
            .padding = qr_size,
            .finders = finders,
            .direction = dir,
            .qr_size = qr_size,
        };
    }
    pub fn deinit(self: FinderPenalty) void {
        self.runEndPositions.deinit();
        self.runHistory.deinit();
    }

    pub fn addHistory(self: *FinderPenalty, current_run_len: usize) !void {
        try self.runHistory.insert(0, current_run_len + self.padding);
        self.padding = 0;
        self.position += current_run_len;
        try self.runEndPositions.insert(0, self.position);
    }
    pub fn countAndAddPatterns(self: *FinderPenalty) !usize {
        const hist = self.runHistory.items;
        const n = hist[1];
        if (n > self.qr_size * 3) return error.WrongSize;

        const core = n > 0 and hist[2] == n and hist[3] == n * 3 and hist[4] == n and hist[5] == n;
        const core_start = self.runEndPositions.items[6];
        const core_end = self.runEndPositions[1];

        var res: usize = 0;
        if (core and hist[0] >= n * 4 and hist[6] >= n) {
            res += 1;
            const start = core_start - n;
            const end = core_end + n * 4;
            try self.finders.append(if (self.direction == .horiz)
                .init(start, self.outer, end - start)
            else
                .init(self.outer, start, end - start));
        }
        return res;
    }
    pub fn terminateAndCount(
        self: *FinderPenalty,
        current_run_color: bool,
        current_run_len: usize,
    ) usize {
        var cur_run = current_run_len;
        if (current_run_color) {
            try self.addHistory(current_run_len);
            cur_run = 0;
        }
        self.padding = self.qr_size;
        try self.addHistory(cur_run);
        if (self.position != self.qr_size) return error.OutOfBounds;
        return try self.countAndAddPatterns();
    }
};
