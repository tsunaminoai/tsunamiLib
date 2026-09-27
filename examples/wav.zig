//! Synthesize a tone, write it as a WAV file into an in-memory buffer,
//! parse it back and verify samples. Also round-trips a struct through
//! io.blob.Blob and shows a corrupted byte being rejected.
const std = @import("std");
const ts = @import("tsunami");
const wav = ts.io.wav;
const blob = ts.io.blob;

const fs = 8000;
const n = 400;

pub fn main(init: std.process.Init) !void {
    var stdout_buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &stdout_buf);
    const w = &stdout.interface;
    const gpa = init.gpa;

    var tone: [n]f32 = undefined;
    for (&tone, 0..) |*s, i| {
        const t = @as(f32, @floatFromInt(i)) / fs;
        s.* = 0.8 * @sin(2 * std.math.pi * 440.0 * t);
    }

    var aw: std.Io.Writer.Allocating = .init(gpa);
    defer aw.deinit();
    const data_size: u32 = n * 2;
    try wav.writeHeader(&aw.writer, .pcm, 1, fs, 16, data_size);
    try wav.writeSamples(f32, &aw.writer, .pcm, 16, &tone);

    var reader: std.Io.Reader = .fixed(aw.written());
    const hdr = try wav.Header.parse(&reader);
    var back: [n]f32 = undefined;
    const got = try wav.readSamples(f32, &reader, hdr, &back);

    var max_err: f32 = 0;
    for (tone, back) |want, have| max_err = @max(max_err, @abs(want - have));

    try w.print("wrote {d} samples @ {d} Hz 16-bit PCM ({d} bytes), read back {d}\n", .{ n, hdr.sample_rate, aw.written().len, got });
    try w.print("max sample error after round-trip: {e:.2}\n", .{max_err});

    // io.blob.Blob: versioned fixed-size struct save/load with a checksum.
    const Settings = struct { sample_rate: u32, gain_pct: u8, channels: u8 };
    const Save = blob.Blob(Settings, .{ 'W', 'A', 'V', 'S' }, 1);
    const settings: Settings = .{ .sample_rate = fs, .gain_pct = 80, .channels = 1 };

    var save_buf: [Save.size]u8 = undefined;
    Save.encode(settings, &save_buf);
    const loaded = try Save.decode(&save_buf);

    var corrupted = save_buf;
    corrupted[6] ^= 0xFF; // flip a payload byte
    const corrupt_result = Save.decode(&corrupted);

    try w.print("blob round-trip: sample_rate={d} gain={d}% channels={d}\n", .{ loaded.sample_rate, loaded.gain_pct, loaded.channels });
    try w.print("corrupted blob decode: {s}\n", .{if (corrupt_result) |_| "accepted (BUG)" else |e| @errorName(e)});
    try w.flush();

    const ok = got == n and max_err < 0.02 and std.meta.eql(loaded, settings) and
        corrupt_result == error.BadChecksum;
    if (!ok) return error.ExampleFailed;
}
