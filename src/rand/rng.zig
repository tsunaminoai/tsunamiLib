/// Rng — a seedable RNG bundle for deterministic, replayable games: the seed
/// it was built from, a `std.Random.DefaultPrng`, and a `random()` accessor.
/// The caller always supplies the seed (never a wall-clock default) so runs
/// stay reproducible for tests and replay.
const std = @import("std");

pub const Rng = struct {
    /// The seed this Rng was constructed (or last reseeded) from. Kept so a run
    /// can be logged and replayed.
    seed: u64,
    prng: std.Random.DefaultPrng,

    pub fn init(seed: u64) Rng {
        return .{ .seed = seed, .prng = std.Random.DefaultPrng.init(seed) };
    }

    /// A `std.Random` valid for immediate use only: it is a fat pointer into
    /// this Rng, so it dangles across a by-value move — call again after one.
    pub fn random(self: *Rng) std.Random {
        return self.prng.random();
    }

    /// Reset the stream to `seed` (replay, or start a fresh deterministic run).
    pub fn reseed(self: *Rng, seed: u64) void {
        self.seed = seed;
        self.prng = std.Random.DefaultPrng.init(seed);
    }
};

// ── Tests ────────────────────────────────────────────────────────────────

const tst = std.testing;

test Rng {
    var a = Rng.init(0xDEAD_BEEF);
    var b = Rng.init(0xDEAD_BEEF);
    for (0..64) |_| {
        try tst.expectEqual(a.random().int(u64), b.random().int(u64));
    }
    try tst.expectEqual(@as(u64, 0xDEAD_BEEF), a.seed);
}

test "different seeds diverge" {
    var a = Rng.init(1);
    var b = Rng.init(2);
    var any_diff = false;
    for (0..16) |_| {
        if (a.random().int(u64) != b.random().int(u64)) any_diff = true;
    }
    try tst.expect(any_diff);
}

test "reseed replays a stream" {
    var r = Rng.init(0x1234);
    const first = r.random().int(u64);
    _ = r.random().int(u64); // advance
    r.reseed(0x1234);
    try tst.expectEqual(first, r.random().int(u64));
    try tst.expectEqual(@as(u64, 0x1234), r.seed);
}

test "random() feeds shuffle-style consumers deterministically" {
    // The intended usage shape: keep Rng as a field, hand random() to a shuffle
    // at the point of use. Two equally seeded Rngs must shuffle identically.
    var a = Rng.init(0x5EED);
    var b = Rng.init(0x5EED);

    var deck_a: [52]u8 = undefined;
    var deck_b: [52]u8 = undefined;
    for (0..52) |i| {
        deck_a[i] = @intCast(i);
        deck_b[i] = @intCast(i);
    }
    a.random().shuffle(u8, &deck_a);
    b.random().shuffle(u8, &deck_b);
    try tst.expectEqualSlices(u8, &deck_a, &deck_b);
}
