const std = @import("std");

/// Proportional-integral loop filter, designed from a normalised loop
/// bandwidth (loop_bw / symbol_or_sample_rate) and damping factor, per the
/// standard 2nd-order PLL design equations. `k0` is the detector/NCO gain
/// (the error signal is divided by it so `advance` returns a correction
/// directly comparable to the tracked quantity).
pub fn LoopFilter(comptime T: type) type {
    return struct {
        alpha: T, // proportional gain
        beta: T, // integral gain
        integrator: T = 0,

        const Self = @This();

        pub fn init(norm_loop_bw: T, damping: T, k0: T) Self {
            const theta = norm_loop_bw / (damping + 1.0 / (4.0 * damping));
            const d = 1.0 + 2.0 * damping * theta + theta * theta;
            return .{
                .alpha = (4.0 * damping * theta) / d / k0,
                .beta = (4.0 * theta * theta) / d / k0,
            };
        }

        /// Feed one error sample, return the loop's correction output.
        pub fn advance(self: *Self, err: T) T {
            self.integrator += self.beta * err;
            return self.alpha * err + self.integrator;
        }

        pub fn reset(self: *Self) void {
            self.integrator = 0;
        }
    };
}

/// Generic phase-locked loop: a numerically-controlled oscillator (phase
/// accumulator) driven by a `LoopFilter`, decoupled from any particular
/// phase detector — the caller computes the error (Costas, Gardner,
/// pilot-tone, whatever) and calls `advance`.
pub fn Pll(comptime T: type) type {
    return struct {
        loop: LoopFilter(T),
        /// Free-running / centre frequency, radians per sample.
        centre_freq: T,
        /// Current NCO frequency, radians per sample.
        freq: T,
        /// Current NCO phase, radians (wrapped to [-pi, pi) by `advance`).
        phase: T = 0,

        const Self = @This();

        pub fn init(norm_loop_bw: T, damping: T, k0: T, centre_freq: T) Self {
            return .{
                .loop = LoopFilter(T).init(norm_loop_bw, damping, k0),
                .centre_freq = centre_freq,
                .freq = centre_freq,
            };
        }

        /// Feed a phase-detector error (radians), advance the NCO by one
        /// sample, and return the new phase.
        pub fn advance(self: *Self, phase_err: T) T {
            const corr = self.loop.advance(phase_err);
            self.freq = self.centre_freq + corr;
            self.phase += self.freq;
            self.phase = wrap(self.phase);
            return self.phase;
        }

        pub fn reset(self: *Self) void {
            self.loop.reset();
            self.freq = self.centre_freq;
            self.phase = 0;
        }

        fn wrap(p: T) T {
            var x = p;
            const two_pi: T = 2.0 * std.math.pi;
            while (x >= std.math.pi) x -= two_pi;
            while (x < -std.math.pi) x += two_pi;
            return x;
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "LoopFilter: step error grows via integration" {
    var lf = LoopFilter(f32).init(0.01, 0.707, 1.0);
    const out0 = lf.advance(0.1);
    try testing.expect(out0 > 0);
    var out: f32 = out0;
    for (0..20) |_| out = lf.advance(0.1);
    try testing.expect(out > out0);
}

test "LoopFilter: zero error gives zero correction after reset" {
    var lf = LoopFilter(f32).init(0.01, 0.707, 1.0);
    _ = lf.advance(0.5);
    lf.reset();
    try testing.expectEqual(@as(f32, 0.0), lf.advance(0.0));
}

test Pll {
    // Detector: error = desired_phase - nco_phase (small-angle, wrapped).
    var pll = Pll(f32).init(0.05, 0.707, 1.0, 0.0);
    const target: f32 = 0.3; // radians, held constant each step
    var last_err: f32 = 1e9;
    for (0..500) |_| {
        const err = target - pll.phase;
        last_err = @abs(err);
        _ = pll.advance(err);
    }
    try testing.expect(last_err < 0.01);
}

test "Pll: reset returns to centre frequency and zero phase" {
    var pll = Pll(f32).init(0.05, 0.707, 1.0, 0.1);
    for (0..50) |_| _ = pll.advance(0.2);
    pll.reset();
    try testing.expectEqual(@as(f32, 0.0), pll.phase);
    try testing.expectEqual(@as(f32, 0.1), pll.freq);
}
