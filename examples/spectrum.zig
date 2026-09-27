//! Two tones → comptime-designed FIR lowpass → windowed FFT; then a biquad
//! highpass on the same input. Shows where each filter puts the energy.
const std = @import("std");
const ts = @import("tsunami");
const dsp = ts.dsp;

const fs = 16_000.0;
const n = 1024;
const F = dsp.fft.Fft(f32, n);

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    var x: [n]f32 = undefined;
    for (&x, 0..) |*v, i| {
        const t = @as(f32, @floatFromInt(i)) / fs;
        v.* = @sin(2 * std.math.pi * 500 * t) + @sin(2 * std.math.pi * 6000 * t);
    }

    // 63 taps, cutoff 2 kHz: computed entirely at compile time into rodata.
    const taps = comptime dsp.fir.lowpass(f32, 63, 2000.0 / fs, .blackman);
    var lp: dsp.fir.Fir(f32, 63) = .init(taps);
    var y_lp: [n]f32 = undefined;
    lp.processBlock(&x, &y_lp);

    var hp = dsp.biquad.Biquad(f32).highpass(2000, fs, std.math.sqrt1_2);
    var y_hp: [n]f32 = undefined;
    for (x, &y_hp) |s, *o| o.* = hp.process(s);

    try w.print("{s:>10} {s:>9} {s:>9}\n", .{ "signal", "500 Hz", "6 kHz" });
    inline for (.{ .{ "input", &x }, .{ "fir lp", &y_lp }, .{ "biquad hp", &y_hp } }) |row| {
        const lo, const hi = tonesDb(row[1]);
        try w.print("{s:>10} {d:>7.1}dB {d:>7.1}dB\n", .{ row[0], lo, hi });
    }
    try w.flush();
}

fn tonesDb(sig: *const [n]f32) struct { f32, f32 } {
    const win = comptime dsp.window.hann(f32, n);
    var s = sig.*;
    dsp.window.apply(f32, &win, &s);
    var scratch: [n]F.C = undefined;
    var spec: [n / 2 + 1]F.C = undefined;
    F.real(&s, &scratch, &spec);
    const bin = struct {
        fn db(sp: []const F.C, hz: f32) f32 {
            const k: usize = @intFromFloat(@round(hz * n / fs));
            return 20 * std.math.log10(sp[k].magnitude() / (n / 4) + 1e-9);
        }
    };
    return .{ bin.db(&spec, 500), bin.db(&spec, 6000) };
}
