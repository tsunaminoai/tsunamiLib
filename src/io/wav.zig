//! RIFF/WAVE reader and writer over `*std.Io.Reader` / `*std.Io.Writer`.
//!
//! Supports PCM 8/16/24/32-bit integer and 32-bit IEEE float samples, any
//! channel count. Untrusted input is never trusted for control flow: chunk
//! sizes are bounds-checked and unknown chunks are skipped rather than
//! rejected, so a malformed file returns an error instead of panicking or
//! reading out of bounds.

const std = @import("std");
const contract = @import("../contract.zig");

/// Sample encoding recognized in the `fmt ` chunk (WAVE_FORMAT_EXTENSIBLE is
/// resolved down to one of these by `Header.parse`).
pub const Format = enum {
    pcm,
    ieee_float,
};

pub const ParseError = error{
    ReadFailed,
    EndOfStream,
    NotRiff,
    NotWave,
    MissingFmtChunk,
    MissingDataChunk,
    ChunkTooSmall,
    InvalidChannelCount,
    InvalidBitsPerSample,
    UnsupportedFormat,
    Overflow,
};

/// Parsed `fmt ` + `data` chunk metadata. Byte order of the underlying stream
/// is always little-endian per the RIFF spec; `readSamples` handles the
/// conversion to host order.
pub const Header = struct {
    format: Format,
    channels: u16,
    sample_rate: u32,
    bits_per_sample: u16,
    /// Size in bytes of the `data` chunk (i.e. all interleaved sample bytes).
    data_size: u32,

    /// Bytes per single-channel sample (e.g. 2 for 16-bit PCM).
    pub fn bytesPerSample(self: Header) u32 {
        return self.bits_per_sample / 8;
    }

    /// Total interleaved samples (channels * frames) in the data chunk.
    pub fn sampleCount(self: Header) u32 {
        const bps = self.bytesPerSample();
        if (bps == 0) return 0;
        return self.data_size / bps;
    }

    /// Frames (per-channel sample groups) in the data chunk.
    pub fn frameCount(self: Header) u32 {
        if (self.channels == 0) return 0;
        return self.sampleCount() / self.channels;
    }

    /// Parse RIFF/WAVE/fmt/data chunk headers from `r`, skipping unknown
    /// chunks. Leaves the reader positioned at the start of sample data.
    /// Never panics or indexes out of bounds on malformed/truncated input.
    pub fn parse(r: *std.Io.Reader) ParseError!Header {
        var riff_id: [4]u8 = undefined;
        try readNoEof(r, &riff_id);
        if (!std.mem.eql(u8, &riff_id, "RIFF")) return error.NotRiff;

        _ = try takeU32(r); // RIFF chunk size (total file size - 8); unused, we stream chunk-by-chunk.

        var wave_id: [4]u8 = undefined;
        try readNoEof(r, &wave_id);
        if (!std.mem.eql(u8, &wave_id, "WAVE")) return error.NotWave;

        var fmt: ?RawFmt = null;
        var data_size: ?u32 = null;

        while (true) {
            var chunk_id: [4]u8 = undefined;
            readNoEof(r, &chunk_id) catch |err| switch (err) {
                error.EndOfStream => break, // no more chunks; require fmt+data below
                else => |e| return e,
            };
            const chunk_size = try takeU32(r);

            if (std.mem.eql(u8, &chunk_id, "fmt ")) {
                fmt = try parseFmtChunk(r, chunk_size);
            } else if (std.mem.eql(u8, &chunk_id, "data")) {
                data_size = chunk_size;
                break; // sample data follows immediately; stop here
            } else {
                try skipChunk(r, chunk_size);
            }
        }

        const f = fmt orelse return error.MissingFmtChunk;
        const size = data_size orelse return error.MissingDataChunk;

        try validateFmt(f);

        return .{
            .format = f.format,
            .channels = f.channels,
            .sample_rate = f.sample_rate,
            .bits_per_sample = f.bits,
            .data_size = size,
        };
    }
};

const RawFmt = struct {
    format: Format,
    channels: u16,
    sample_rate: u32,
    bits: u16,
};

fn parseFmtChunk(r: *std.Io.Reader, chunk_size: u32) ParseError!RawFmt {
    // Minimal PCM fmt chunk: code(2) channels(2) rate(4) byte_rate(4) align(2) bits(2) = 16
    if (chunk_size < 16) return error.ChunkTooSmall;

    const code = try takeU16(r);
    const channels = try takeU16(r);
    const sample_rate = try takeU32(r);
    _ = try takeU32(r); // byte rate; derivable, not load-bearing here
    _ = try takeU16(r); // block align; derivable, not load-bearing here
    const bits = try takeU16(r);

    var remaining: u32 = chunk_size - 16;
    var resolved_code = code;

    if (code == 0xFFFE) {
        // WAVE_FORMAT_EXTENSIBLE: cbSize(2) valid_bits(2) channel_mask(4)
        // then a 16-byte SubFormat GUID whose first 2 bytes are the real
        // format code.
        if (remaining < 24) return error.ChunkTooSmall;
        _ = try takeU16(r); // cbSize
        _ = try takeU16(r); // valid bits per sample
        _ = try takeU32(r); // channel mask
        resolved_code = try takeU16(r); // SubFormat GUID leading u16
        try r.discardAll(14); // rest of the GUID
        remaining -= 24;
    }

    if (remaining > 0) try r.discardAll(remaining);

    const format: Format = switch (resolved_code) {
        1 => .pcm,
        3 => .ieee_float,
        else => return error.UnsupportedFormat,
    };

    return .{ .format = format, .channels = channels, .sample_rate = sample_rate, .bits = bits };
}

fn validateFmt(f: RawFmt) ParseError!void {
    if (f.channels == 0) return error.InvalidChannelCount;
    switch (f.bits) {
        8, 16, 24, 32 => {},
        else => return error.InvalidBitsPerSample,
    }
    if (f.format == .ieee_float and f.bits != 32) return error.UnsupportedFormat;
}

fn skipChunk(r: *std.Io.Reader, chunk_size: u32) ParseError!void {
    // RIFF chunks are word-aligned: an odd-sized chunk carries one pad byte
    // not counted in its declared size.
    const padded = chunk_size + (chunk_size & 1);
    try r.discardAll(padded);
}

fn readNoEof(r: *std.Io.Reader, buf: []u8) ParseError!void {
    r.readSliceAll(buf) catch |err| switch (err) {
        error.ReadFailed => return error.ReadFailed,
        error.EndOfStream => return error.EndOfStream,
    };
}

fn takeU32(r: *std.Io.Reader) ParseError!u32 {
    return r.takeInt(u32, .little) catch |err| switch (err) {
        error.ReadFailed => return error.ReadFailed,
        error.EndOfStream => return error.EndOfStream,
    };
}

fn takeU16(r: *std.Io.Reader) ParseError!u16 {
    return r.takeInt(u16, .little) catch |err| switch (err) {
        error.ReadFailed => return error.ReadFailed,
        error.EndOfStream => return error.EndOfStream,
    };
}

pub const ReadSamplesError = std.Io.Reader.Error || error{EndOfStream};

/// Read and normalize up to `hdr.sampleCount()` interleaved samples from `r`
/// into `out`, converting integer PCM or float32 source samples to `T`
/// (typically `f32`, normalized to `[-1, 1]`). Multi-channel samples stay
/// interleaved: frame `t`, channel `c` lands at `out[t * channels + c]`.
/// Returns the number of samples actually read (may be less than
/// `hdr.sampleCount()` if `out` is smaller, or if the stream is short).
pub fn readSamples(comptime T: type, r: *std.Io.Reader, hdr: Header, out: []T) ReadSamplesError!usize {
    contract.require(out.len > 0 or hdr.sampleCount() == 0, "wav.readSamples: out.len == 0 with samples pending");

    const want = @min(out.len, hdr.sampleCount());
    var i: usize = 0;
    switch (hdr.format) {
        .pcm => switch (hdr.bits_per_sample) {
            8 => while (i < want) : (i += 1) {
                out[i] = normalizeU8(try r.takeByte());
            },
            16 => while (i < want) : (i += 1) {
                out[i] = normalizeInt(T, i16, try r.takeInt(i16, .little));
            },
            24 => while (i < want) : (i += 1) {
                out[i] = normalizeInt(T, i24, try readI24(r));
            },
            32 => while (i < want) : (i += 1) {
                out[i] = normalizeInt(T, i32, try r.takeInt(i32, .little));
            },
            else => unreachable, // validated in Header.parse
        },
        .ieee_float => while (i < want) : (i += 1) {
            const bits = try r.takeInt(u32, .little);
            out[i] = @floatCast(@as(f32, @bitCast(bits)));
        },
    }
    return i;
}

fn readI24(r: *std.Io.Reader) std.Io.Reader.Error!i24 {
    const b = try r.takeArray(3);
    const u: u24 = @as(u24, b[0]) | (@as(u24, b[1]) << 8) | (@as(u24, b[2]) << 16);
    return @bitCast(u);
}

fn normalizeU8(v: u8) f32 {
    // 8-bit PCM is conventionally unsigned with 128 as the zero point.
    return (@as(f32, @floatFromInt(v)) - 128.0) / 128.0;
}

fn normalizeInt(comptime T: type, comptime S: type, v: S) T {
    const max_mag: f32 = @floatFromInt(-@as(i64, std.math.minInt(S)));
    const f: f32 = @as(f32, @floatFromInt(v)) / max_mag;
    return @floatCast(f);
}

/// Write a canonical 44-byte PCM/IEEE-float WAV header (single `fmt ` chunk,
/// no extensible extension) via `w`. `bits_per_sample` must be one of
/// 8/16/24/32; `format` must be `.ieee_float` only when `bits_per_sample`
/// is 32. `data_size` is the total byte length of the sample data that will
/// follow (see `writeSamples`).
pub fn writeHeader(
    w: *std.Io.Writer,
    format: Format,
    channels: u16,
    sample_rate: u32,
    bits_per_sample: u16,
    data_size: u32,
) std.Io.Writer.Error!void {
    contract.require(channels > 0, "wav.writeHeader: channels == 0");
    contract.require(
        bits_per_sample == 8 or bits_per_sample == 16 or bits_per_sample == 24 or bits_per_sample == 32,
        "wav.writeHeader: bits_per_sample must be 8/16/24/32",
    );
    contract.require(
        format != .ieee_float or bits_per_sample == 32,
        "wav.writeHeader: ieee_float requires bits_per_sample == 32",
    );

    const block_align: u16 = channels * (bits_per_sample / 8);
    const byte_rate: u32 = sample_rate * @as(u32, block_align);
    const fmt_code: u16 = switch (format) {
        .pcm => 1,
        .ieee_float => 3,
    };
    const riff_size: u32 = 36 + data_size; // 4("WAVE") + (8+16) fmt + (8) data hdr + data

    try w.writeAll("RIFF");
    try w.writeInt(u32, riff_size, .little);
    try w.writeAll("WAVE");

    try w.writeAll("fmt ");
    try w.writeInt(u32, 16, .little);
    try w.writeInt(u16, fmt_code, .little);
    try w.writeInt(u16, channels, .little);
    try w.writeInt(u32, sample_rate, .little);
    try w.writeInt(u32, byte_rate, .little);
    try w.writeInt(u16, block_align, .little);
    try w.writeInt(u16, bits_per_sample, .little);

    try w.writeAll("data");
    try w.writeInt(u32, data_size, .little);
}

/// Write `samples` (interleaved, values expected in `[-1, 1]` for PCM
/// formats) to `w`, converting from `T` to `bits_per_sample`-wide PCM or to
/// float32, matching a header written by `writeHeader`. Out-of-range PCM
/// values are clamped rather than wrapped.
pub fn writeSamples(
    comptime T: type,
    w: *std.Io.Writer,
    format: Format,
    bits_per_sample: u16,
    samples: []const T,
) std.Io.Writer.Error!void {
    contract.require(
        bits_per_sample == 8 or bits_per_sample == 16 or bits_per_sample == 24 or bits_per_sample == 32,
        "wav.writeSamples: bits_per_sample must be 8/16/24/32",
    );

    switch (format) {
        .pcm => switch (bits_per_sample) {
            8 => for (samples) |s| try w.writeByte(denormalizeU8(s)),
            16 => for (samples) |s| try w.writeInt(i16, denormalizeInt(i16, s), .little),
            24 => for (samples) |s| try writeI24(w, denormalizeInt(i24, s)),
            32 => for (samples) |s| try w.writeInt(i32, denormalizeInt(i32, s), .little),
            else => unreachable, // checked above
        },
        .ieee_float => for (samples) |s| {
            const f: f32 = @floatCast(s);
            try w.writeInt(u32, @bitCast(f), .little);
        },
    }
}

fn writeI24(w: *std.Io.Writer, v: i24) std.Io.Writer.Error!void {
    const u: u24 = @bitCast(v);
    try w.writeByte(@truncate(u));
    try w.writeByte(@truncate(u >> 8));
    try w.writeByte(@truncate(u >> 16));
}

fn denormalizeU8(v: anytype) u8 {
    const f: f32 = @floatCast(std.math.clamp(v, -1.0, 1.0));
    const scaled = (f + 1.0) * 128.0;
    return @intFromFloat(std.math.clamp(scaled, 0.0, 255.0));
}

fn denormalizeInt(comptime S: type, v: anytype) S {
    const f: f32 = @floatCast(std.math.clamp(v, -1.0, 1.0));
    const max_mag: f32 = @floatFromInt(-@as(i64, std.math.minInt(S)));
    const scaled = f * max_mag;
    const lo: f32 = @floatFromInt(std.math.minInt(S));
    const hi: f32 = @floatFromInt(std.math.maxInt(S));
    const clamped = std.math.clamp(scaled, lo, hi);
    return @intFromFloat(clamped);
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "fuzz parse" {
    try std.testing.fuzz({}, struct {
        fn f(_: void, s: *std.testing.Smith) !void {
            var buf: [512]u8 = undefined;
            const n = s.slice(&buf);
            var r: std.Io.Reader = .fixed(buf[0..n]);
            _ = Header.parse(&r) catch return;
        }
    }.f, .{});
}

test "parse rejects non-RIFF" {
    var r: std.Io.Reader = .fixed("JUNKxxxxWAVE");
    try testing.expectError(error.NotRiff, Header.parse(&r));
}

test "parse rejects non-WAVE" {
    var buf: [12]u8 = undefined;
    @memcpy(buf[0..4], "RIFF");
    std.mem.writeInt(u32, buf[4..8], 4, .little);
    @memcpy(buf[8..12], "JUNK");
    var r: std.Io.Reader = .fixed(&buf);
    try testing.expectError(error.NotWave, Header.parse(&r));
}

test "parse skips unknown chunks before fmt/data" {
    var out_buf: [256]u8 = undefined;
    var w: std.Io.Writer = .fixed(&out_buf);
    try w.writeAll("RIFF");
    try w.writeInt(u32, 0, .little); // riff size, unused by parse
    try w.writeAll("WAVE");

    // Unknown "JUNK" chunk with odd size 3 (+1 pad byte).
    try w.writeAll("JUNK");
    try w.writeInt(u32, 3, .little);
    try w.writeAll("abc");
    try w.writeByte(0); // pad

    // fmt chunk (16-bit mono PCM @ 8000 Hz).
    try w.writeAll("fmt ");
    try w.writeInt(u32, 16, .little);
    try w.writeInt(u16, 1, .little); // PCM
    try w.writeInt(u16, 1, .little); // channels
    try w.writeInt(u32, 8000, .little); // sample rate
    try w.writeInt(u32, 8000 * 2, .little); // byte rate
    try w.writeInt(u16, 2, .little); // block align
    try w.writeInt(u16, 16, .little); // bits per sample

    // data chunk.
    try w.writeAll("data");
    try w.writeInt(u32, 4, .little);
    try writeSamples(f32, &w, .pcm, 16, &.{ 0.5, -0.5 });

    var r: std.Io.Reader = .fixed(w.buffered());
    const hdr = try Header.parse(&r);
    try testing.expectEqual(Format.pcm, hdr.format);
    try testing.expectEqual(@as(u16, 1), hdr.channels);
    try testing.expectEqual(@as(u32, 8000), hdr.sample_rate);
    try testing.expectEqual(@as(u16, 16), hdr.bits_per_sample);
    try testing.expectEqual(@as(u32, 4), hdr.data_size);
}

test "readSamples requires nonzero out when samples pending" {
    // contract.require would panic; instead check the escape hatch: zero
    // samples pending is fine with an empty out slice.
    var buf: [64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    try writeHeader(&w, .pcm, 1, 8000, 16, 0);
    var r: std.Io.Reader = .fixed(w.buffered());
    const hdr = try Header.parse(&r);
    var out: [0]f32 = undefined;
    const n = try readSamples(f32, &r, hdr, &out);
    try testing.expectEqual(@as(usize, 0), n);
}

test Header {
    var buf: [256]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    const samples = [_]f32{ 0.0, 0.5, -0.5, 1.0, -1.0 };
    const data_size: u32 = samples.len * 2;
    try writeHeader(&w, .pcm, 1, 44100, 16, data_size);
    try writeSamples(f32, &w, .pcm, 16, &samples);

    var r: std.Io.Reader = .fixed(w.buffered());
    const hdr = try Header.parse(&r);
    try testing.expectEqual(Format.pcm, hdr.format);
    try testing.expectEqual(@as(u16, 1), hdr.channels);
    try testing.expectEqual(@as(u32, 44100), hdr.sample_rate);
    try testing.expectEqual(@as(u16, 16), hdr.bits_per_sample);
    try testing.expectEqual(@as(u32, data_size), hdr.data_size);
    try testing.expectEqual(@as(u32, samples.len), hdr.sampleCount());

    var out: [samples.len]f32 = undefined;
    const n = try readSamples(f32, &r, hdr, &out);
    try testing.expectEqual(samples.len, n);
    for (samples, out) |want, got| try testing.expectApproxEqAbs(want, got, 1e-3);
}

test "round trip: stereo 24-bit PCM" {
    var buf: [256]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    // Interleaved L/R frames.
    const samples = [_]f32{ 0.25, -0.25, 0.75, -0.75 };
    const data_size: u32 = samples.len * 3;
    try writeHeader(&w, .pcm, 2, 48000, 24, data_size);
    try writeSamples(f32, &w, .pcm, 24, &samples);

    var r: std.Io.Reader = .fixed(w.buffered());
    const hdr = try Header.parse(&r);
    try testing.expectEqual(@as(u16, 2), hdr.channels);
    try testing.expectEqual(@as(u16, 24), hdr.bits_per_sample);
    try testing.expectEqual(@as(u32, 2), hdr.frameCount());

    var out: [samples.len]f32 = undefined;
    const n = try readSamples(f32, &r, hdr, &out);
    try testing.expectEqual(samples.len, n);
    for (samples, out) |want, got| try testing.expectApproxEqAbs(want, got, 1e-2);
}

test "round trip: float32" {
    var buf: [256]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    const samples = [_]f32{ 0.123, -0.456, 0.999 };
    const data_size: u32 = samples.len * 4;
    try writeHeader(&w, .ieee_float, 1, 96000, 32, data_size);
    try writeSamples(f32, &w, .ieee_float, 32, &samples);

    var r: std.Io.Reader = .fixed(w.buffered());
    const hdr = try Header.parse(&r);
    try testing.expectEqual(Format.ieee_float, hdr.format);

    var out: [samples.len]f32 = undefined;
    const n = try readSamples(f32, &r, hdr, &out);
    try testing.expectEqual(samples.len, n);
    for (samples, out) |want, got| try testing.expectApproxEqAbs(want, got, 1e-6);
}

test "round trip: 8-bit PCM" {
    var buf: [64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    const samples = [_]f32{ 0.0, 1.0, -1.0 };
    try writeHeader(&w, .pcm, 1, 8000, 8, samples.len);
    try writeSamples(f32, &w, .pcm, 8, &samples);

    var r: std.Io.Reader = .fixed(w.buffered());
    const hdr = try Header.parse(&r);
    try testing.expectEqual(@as(u16, 8), hdr.bits_per_sample);

    var out: [samples.len]f32 = undefined;
    const n = try readSamples(f32, &r, hdr, &out);
    try testing.expectEqual(samples.len, n);
    try testing.expectApproxEqAbs(@as(f32, 0.0), out[0], 0.02);
    try testing.expectApproxEqAbs(@as(f32, 1.0), out[1], 0.02);
    try testing.expectApproxEqAbs(@as(f32, -1.0), out[2], 0.02);
}

test "readSamples truncates to out.len" {
    var buf: [256]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);
    const samples = [_]f32{ 0.1, 0.2, 0.3, 0.4 };
    try writeHeader(&w, .pcm, 1, 8000, 16, samples.len * 2);
    try writeSamples(f32, &w, .pcm, 16, &samples);

    var r: std.Io.Reader = .fixed(w.buffered());
    const hdr = try Header.parse(&r);
    var out: [2]f32 = undefined;
    const n = try readSamples(f32, &r, hdr, &out);
    try testing.expectEqual(@as(usize, 2), n);
}

// Fuzz mode is broken in the 0.16.0 test runner; mutate a valid file instead.
test "mutated files never crash parse or readSamples" {
    var file: [44 + 64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&file);
    try writeHeader(&w, .pcm, 2, 8000, 16, 64);
    const src = [_]f32{ 0.5, -0.5 } ** 16;
    try writeSamples(f32, &w, .pcm, 16, &src);
    var prng: std.Random.DefaultPrng = .init(0x3a7);
    const r = prng.random();
    var out: [64]f32 = undefined;
    for (0..20_000) |_| {
        var m = file;
        for (0..r.intRangeAtMost(usize, 1, 6)) |_| m[r.uintLessThan(usize, m.len)] = r.int(u8);
        const len = r.uintAtMost(usize, m.len);
        var rd: std.Io.Reader = .fixed(m[0..len]);
        const hdr = Header.parse(&rd) catch continue;
        _ = readSamples(f32, &rd, hdr, &out) catch continue;
    }
}
