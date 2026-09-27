//! Random bytes → 64-QAM pack/map → AWGN → LLRs → LDPC(1944,r5/6) decode.
//! Shows pre-FEC bit errors from the noisy channel collapsing to zero post-FEC.
const std = @import("std");
const ts = @import("tsunami");
const coding = ts.coding;
const dsp = ts.dsp;

const Q = coding.constellation.Qam(f32, 6);
const C = coding.ldpc.Ldpc1944R56;
const snr_db: f32 = 18.0;

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    var prng: std.Random.DefaultPrng = .init(0xC0FFEE);
    var rand = prng.random();

    var info: [C.k]u1 = undefined;
    for (&info) |*b| b.* = rand.int(u1);
    var cw: [C.n]u1 = undefined;
    C.encode(&info, &cw);

    // Pack codeword bits (n a multiple of 6) into 6-bit QAM symbols.
    const n_syms = C.n / 6;
    var syms: [n_syms]Q.C = undefined;
    for (0..n_syms) |i| {
        var code: Q.Code = 0;
        for (0..6) |b| code = (code << 1) | cw[i * 6 + b];
        syms[i] = Q.map(code);
    }

    // AWGN at the target Es/N0; noise_var derived from the same sigma the
    // channel actually added, so the LLRs are calibrated (not just signed).
    var re: [n_syms]f32 = undefined;
    var im: [n_syms]f32 = undefined;
    for (syms, &re, &im) |s, *r, *i| {
        r.* = s.re;
        i.* = s.im;
    }
    dsp.channel.awgn(&re, snr_db, &rand);
    dsp.channel.awgn(&im, snr_db, &rand);
    const noise_var = std.math.pow(f32, 10.0, -snr_db / 10.0);

    var llr: [C.n]f32 = undefined;
    var pre_fec_errors: usize = 0;
    for (0..n_syms) |i| {
        const s: Q.C = .{ .re = re[i], .im = im[i] };
        var out: [6]f32 = undefined;
        Q.llr(s, noise_var, &out);
        @memcpy(llr[i * 6 ..][0..6], &out);
        const hard = Q.slice(s);
        var want: Q.Code = 0;
        for (0..6) |b| want = (want << 1) | cw[i * 6 + b];
        pre_fec_errors += @popCount(hard ^ want);
    }

    var decoded: [C.n]u1 = undefined;
    const ok = C.decode(&llr, &decoded, coding.ldpc.default_iters);
    var post_fec_errors: usize = 0;
    for (decoded, cw) |d, c| post_fec_errors += @intFromBool(d != c);

    try w.print("QAM-64 + LDPC(1944,r5/6) over AWGN @ {d:.1} dB Es/N0\n", .{snr_db});
    try w.print("pre-FEC bit errors:  {d} / {d}\n", .{ pre_fec_errors, C.n });
    try w.print("post-FEC bit errors: {d} / {d}  (decode {s})\n", .{
        post_fec_errors, C.n, if (ok) "converged" else "failed",
    });
    try w.flush();

    if (!ok or post_fec_errors != 0) return error.ExampleFailed;
}
