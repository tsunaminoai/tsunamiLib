const std = @import("std");

/// RBJ cookbook biquad, transposed direct form II (2 state vars, numerically
/// robust). Coefficients are pre-normalised so `a0 == 1` is implicit.
pub fn Biquad(comptime T: type) type {
    return struct {
        b0: T,
        b1: T,
        b2: T,
        a1: T,
        a2: T,
        s1: T = 0,
        s2: T = 0,

        const Self = @This();

        fn make(b0: T, b1: T, b2: T, a0: T, a1: T, a2: T) Self {
            return .{
                .b0 = b0 / a0,
                .b1 = b1 / a0,
                .b2 = b2 / a0,
                .a1 = a1 / a0,
                .a2 = a2 / a0,
            };
        }

        /// fc: cutoff Hz, fs: sample rate Hz, q: quality factor (0.7071 ≈ Butterworth).
        pub fn lowpass(fc: T, fs: T, q: T) Self {
            const c = Coeffs(T).init(fc, fs, q);
            const b1 = 1.0 - c.cs;
            return make(b1 / 2, b1, b1 / 2, 1.0 + c.alpha, -2.0 * c.cs, 1.0 - c.alpha);
        }

        pub fn highpass(fc: T, fs: T, q: T) Self {
            const c = Coeffs(T).init(fc, fs, q);
            const b1 = -(1.0 + c.cs);
            return make(-b1 / 2, b1, -b1 / 2, 1.0 + c.alpha, -2.0 * c.cs, 1.0 - c.alpha);
        }

        pub fn bandpass(fc: T, fs: T, q: T) Self {
            const c = Coeffs(T).init(fc, fs, q);
            return make(c.alpha, 0, -c.alpha, 1.0 + c.alpha, -2.0 * c.cs, 1.0 - c.alpha);
        }

        pub fn notch(fc: T, fs: T, q: T) Self {
            const c = Coeffs(T).init(fc, fs, q);
            return make(1.0, -2.0 * c.cs, 1.0, 1.0 + c.alpha, -2.0 * c.cs, 1.0 - c.alpha);
        }

        /// db_gain: boost/cut at fc, in dB.
        pub fn peak(fc: T, fs: T, q: T, db_gain: T) Self {
            const c = Coeffs(T).init(fc, fs, q);
            const a = std.math.pow(T, 10.0, db_gain / 40.0);
            return make(1.0 + c.alpha * a, -2.0 * c.cs, 1.0 - c.alpha * a, 1.0 + c.alpha / a, -2.0 * c.cs, 1.0 - c.alpha / a);
        }

        pub fn lowshelf(fc: T, fs: T, q: T, db_gain: T) Self {
            const c = Coeffs(T).init(fc, fs, q);
            const a = std.math.pow(T, 10.0, db_gain / 40.0);
            const beta = 2.0 * @sqrt(a) * c.alpha;
            return make(
                a * ((a + 1.0) - (a - 1.0) * c.cs + beta),
                2.0 * a * ((a - 1.0) - (a + 1.0) * c.cs),
                a * ((a + 1.0) - (a - 1.0) * c.cs - beta),
                (a + 1.0) + (a - 1.0) * c.cs + beta,
                -2.0 * ((a - 1.0) + (a + 1.0) * c.cs),
                (a + 1.0) + (a - 1.0) * c.cs - beta,
            );
        }

        pub fn highshelf(fc: T, fs: T, q: T, db_gain: T) Self {
            const c = Coeffs(T).init(fc, fs, q);
            const a = std.math.pow(T, 10.0, db_gain / 40.0);
            const beta = 2.0 * @sqrt(a) * c.alpha;
            return make(
                a * ((a + 1.0) + (a - 1.0) * c.cs + beta),
                -2.0 * a * ((a - 1.0) + (a + 1.0) * c.cs),
                a * ((a + 1.0) + (a - 1.0) * c.cs - beta),
                (a + 1.0) - (a - 1.0) * c.cs + beta,
                2.0 * ((a - 1.0) - (a + 1.0) * c.cs),
                (a + 1.0) - (a - 1.0) * c.cs - beta,
            );
        }

        /// Build directly from normalised coefficients (a0 implicit = 1).
        pub fn fromCoeffs(b0: T, b1: T, b2: T, a1: T, a2: T) Self {
            return .{ .b0 = b0, .b1 = b1, .b2 = b2, .a1 = a1, .a2 = a2 };
        }

        /// Transposed direct form II: one multiply-add per state var, no
        /// separate input/output history — fewer roundings than DF-I.
        pub fn process(self: *Self, x: T) T {
            const y = self.b0 * x + self.s1;
            self.s1 = self.b1 * x - self.a1 * y + self.s2;
            self.s2 = self.b2 * x - self.a2 * y;
            return y;
        }

        pub fn reset(self: *Self) void {
            self.s1 = 0;
            self.s2 = 0;
        }
    };
}

fn Coeffs(comptime T: type) type {
    return struct {
        cs: T,
        alpha: T,

        fn init(fc: T, fs: T, q: T) @This() {
            const omega = 2.0 * std.math.pi * fc / fs;
            const sn = @sin(omega);
            const cs = @cos(omega);
            return .{ .cs = cs, .alpha = sn / (2.0 * q) };
        }
    };
}

/// Cascade of `sections` biquads (each section ~12 dB/octave), applied in series.
pub fn Cascade(comptime T: type, comptime sections: usize) type {
    return struct {
        stages: [sections]Biquad(T),

        const Self = @This();

        pub fn init(stages: [sections]Biquad(T)) Self {
            return .{ .stages = stages };
        }

        pub fn process(self: *Self, x: T) T {
            var y = x;
            for (&self.stages) |*s| y = s.process(y);
            return y;
        }

        pub fn reset(self: *Self) void {
            for (&self.stages) |*s| s.reset();
        }
    };
}

/// Butterworth lowpass of the given even or odd `order`, built as a cascade
/// of ⌈order/2⌉ biquad sections. Each section's Q comes from the analog
/// Butterworth pole angles θ_k = π(2k+1)/(2·order), Q_k = 1/(2·cos θ_k)
/// (k = 0..⌈order/2⌉), all sections sharing the same cutoff — the standard
/// even-order factorisation. Odd order gets one extra first-order section
/// folded into a biquad with a2 = 0 (Q → ∞ treated as a first-order RC).
pub fn butterworth(comptime T: type, comptime order: usize, fs: T, fc: T) Cascade(T, (order + 1) / 2) {
    const sections = (order + 1) / 2;
    var stages: [sections]Biquad(T) = undefined;
    const pairs = order / 2;
    for (0..pairs) |k| {
        const theta = std.math.pi * (2.0 * @as(T, @floatFromInt(k)) + 1.0) / (2.0 * @as(T, @floatFromInt(order)));
        const q = 1.0 / (2.0 * @cos(theta));
        stages[k] = Biquad(T).lowpass(fc, fs, q);
    }
    if (order % 2 == 1) {
        // First-order RC lowpass via bilinear transform, expressed as a
        // degenerate biquad (b2 = a2 = 0).
        const k = std.math.tan(std.math.pi * fc / fs);
        const b0 = k / (1.0 + k);
        const a1 = (k - 1.0) / (1.0 + k);
        stages[pairs] = Biquad(T).fromCoeffs(b0, b0, 0, a1, 0);
    }
    return Cascade(T, sections).init(stages);
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

fn magnitudeAt(comptime T: type, filt: anytype, fc_hz: T, fs: T, n: usize) T {
    var f = filt;
    var peak_out: T = 0;
    var settle: usize = 0;
    while (settle < n) : (settle += 1) {
        const t: T = @floatFromInt(settle);
        const x = @sin(2.0 * std.math.pi * fc_hz * t / fs);
        const y = f.process(x);
        if (settle > n / 2) peak_out = @max(peak_out, @abs(y));
    }
    return peak_out;
}

test "lowpass: passes DC-ish low tone, attenuates high tone" {
    var lp = Biquad(f64).lowpass(1000.0, 48000.0, std.math.sqrt1_2);
    const lo = magnitudeAt(f64, lp, 100.0, 48000.0, 4000);
    lp.reset();
    const hi = magnitudeAt(f64, lp, 15000.0, 48000.0, 4000);
    try testing.expect(lo > 0.9);
    try testing.expect(hi < 0.1);
}

test "highpass: attenuates low tone, passes high tone" {
    var hp = Biquad(f64).highpass(1000.0, 48000.0, std.math.sqrt1_2);
    const lo = magnitudeAt(f64, hp, 50.0, 48000.0, 4000);
    hp.reset();
    const hi = magnitudeAt(f64, hp, 15000.0, 48000.0, 4000);
    try testing.expect(lo < 0.1);
    try testing.expect(hi > 0.9);
}

test "bandpass: passes centre, attenuates far away" {
    var bp = Biquad(f64).bandpass(1000.0, 48000.0, 4.0);
    const centre = magnitudeAt(f64, bp, 1000.0, 48000.0, 4000);
    bp.reset();
    const far = magnitudeAt(f64, bp, 8000.0, 48000.0, 4000);
    try testing.expect(centre > 0.5);
    try testing.expect(far < 0.1);
}

test "notch: rejects centre, passes far away" {
    var nf = Biquad(f64).notch(1000.0, 48000.0, 4.0);
    const centre = magnitudeAt(f64, nf, 1000.0, 48000.0, 4000);
    nf.reset();
    const far = magnitudeAt(f64, nf, 8000.0, 48000.0, 4000);
    try testing.expect(centre < 0.1);
    try testing.expect(far > 0.85);
}

test "peak: boosts centre relative to unity gain elsewhere" {
    var pk = Biquad(f64).peak(1000.0, 48000.0, 2.0, 12.0);
    const centre = magnitudeAt(f64, pk, 1000.0, 48000.0, 4000);
    pk.reset();
    const far = magnitudeAt(f64, pk, 100.0, 48000.0, 4000);
    try testing.expect(centre > far);
}

test "lowshelf/highshelf: boost applies at the expected end" {
    const ls = Biquad(f64).lowshelf(1000.0, 48000.0, std.math.sqrt1_2, 12.0);
    const lo = magnitudeAt(f64, ls, 50.0, 48000.0, 4000);
    try testing.expect(lo > 1.5); // ~4x boost

    const hs = Biquad(f64).highshelf(1000.0, 48000.0, std.math.sqrt1_2, 12.0);
    const hi = magnitudeAt(f64, hs, 15000.0, 48000.0, 4000);
    try testing.expect(hi > 1.5);
}

test "Cascade: two lowpass sections attenuate more steeply than one" {
    const one = Biquad(f64).lowpass(1000.0, 48000.0, std.math.sqrt1_2);
    const single = magnitudeAt(f64, one, 4000.0, 48000.0, 4000);

    var casc = Cascade(f64, 2).init(.{ one, one });
    const double = magnitudeAt(f64, casc, 4000.0, 48000.0, 4000);
    _ = &casc;
    try testing.expect(double < single);
}

test "butterworth: order 4 lowpass passes low, attenuates high" {
    var bw = butterworth(f64, 4, 48000.0, 2000.0);
    const lo = magnitudeAt(f64, bw, 200.0, 48000.0, 4000);
    bw.reset();
    const hi = magnitudeAt(f64, bw, 16000.0, 48000.0, 4000);
    try testing.expect(lo > 0.85);
    try testing.expect(hi < 0.02);
}

test "butterworth: odd order (3) also builds and filters" {
    var bw = butterworth(f64, 3, 48000.0, 2000.0);
    const lo = magnitudeAt(f64, bw, 200.0, 48000.0, 4000);
    bw.reset();
    const hi = magnitudeAt(f64, bw, 16000.0, 48000.0, 4000);
    try testing.expect(lo > 0.85);
    try testing.expect(hi < 0.1);
}

test "butterworth: -3dB near cutoff for order 2" {
    const bw = butterworth(f64, 2, 48000.0, 4000.0);
    const at_fc = magnitudeAt(f64, bw, 4000.0, 48000.0, 8000);
    // Expect roughly 0.707 (±0.1) amplitude at cutoff.
    try testing.expect(at_fc > 0.55 and at_fc < 0.85);
}
