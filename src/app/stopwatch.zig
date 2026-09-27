const std = @import("std");
const Io = std.Io;

/// Monotonic (`.awake`) stopwatch with pause and lap support.
pub const Stopwatch = struct {
    start_ts: Io.Timestamp,
    lap_ts: Io.Timestamp,
    paused_at: ?Io.Timestamp = null,

    pub fn start(io: Io) Stopwatch {
        const now = Io.Timestamp.now(io, .awake);
        return .{ .start_ts = now, .lap_ts = now };
    }

    pub fn elapsed(sw: *const Stopwatch, io: Io) Io.Duration {
        const end = sw.paused_at orelse Io.Timestamp.now(io, .awake);
        return sw.start_ts.durationTo(end);
    }

    /// Time since the previous lap (or start), then begins a new lap.
    pub fn lap(sw: *Stopwatch, io: Io) Io.Duration {
        const now = sw.paused_at orelse Io.Timestamp.now(io, .awake);
        defer sw.lap_ts = now;
        return sw.lap_ts.durationTo(now);
    }

    pub fn pause(sw: *Stopwatch, io: Io) void {
        if (sw.paused_at == null) sw.paused_at = Io.Timestamp.now(io, .awake);
    }

    pub fn unpause(sw: *Stopwatch, io: Io) void {
        const p = sw.paused_at orelse return;
        const gap = p.durationTo(Io.Timestamp.now(io, .awake)).nanoseconds;
        sw.start_ts.nanoseconds += gap;
        sw.lap_ts.nanoseconds += gap;
        sw.paused_at = null;
    }
};

/// Pure dt-driven countdown: deterministic, no clock access, frame-loop friendly.
pub fn Countdown(comptime T: type) type {
    return struct {
        const Self = @This();
        duration: T,
        left: T,
        repeat: bool = false,

        pub fn init(duration: T, repeat: bool) Self {
            return .{ .duration = duration, .left = duration, .repeat = repeat };
        }

        /// True on the tick it fires. Repeating timers carry overshoot forward.
        pub fn update(c: *Self, dt: T) bool {
            if (c.left <= 0 and !c.repeat) return false;
            c.left -= dt;
            if (c.left > 0) return false;
            if (c.repeat) c.left += c.duration;
            return true;
        }

        pub fn reset(c: *Self) void {
            c.left = c.duration;
        }

        pub fn progress(c: *const Self) T {
            return std.math.clamp(1 - c.left / c.duration, 0, 1);
        }
    };
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test Stopwatch {
    const io = testing.io;
    var sw: Stopwatch = .start(io);
    try io.sleep(.fromMilliseconds(2), .awake);
    try testing.expect(sw.elapsed(io).toMilliseconds() >= 2);
    const l1 = sw.lap(io);
    try testing.expect(l1.toMilliseconds() >= 2);
    try testing.expect(sw.lap(io).nanoseconds < l1.nanoseconds);

    sw.pause(io);
    const frozen = sw.elapsed(io);
    try io.sleep(.fromMilliseconds(2), .awake);
    try testing.expectEqual(frozen.nanoseconds, sw.elapsed(io).nanoseconds);
    sw.unpause(io);
    try testing.expect(sw.elapsed(io).nanoseconds >= frozen.nanoseconds);
    try testing.expect(sw.elapsed(io).nanoseconds < frozen.nanoseconds + std.time.ns_per_ms);
}

test "countdown one-shot and repeating" {
    var c: Countdown(f32) = .init(1.0, false);
    try testing.expect(!c.update(0.6));
    try testing.expectApproxEqAbs(@as(f32, 0.6), c.progress(), 1e-6);
    try testing.expect(c.update(0.6));
    try testing.expect(!c.update(0.6));

    var r: Countdown(f32) = .init(1.0, true);
    var fires: usize = 0;
    for (0..10) |_| fires += @intFromBool(r.update(0.25));
    try testing.expectEqual(@as(usize, 2), fires);
}
