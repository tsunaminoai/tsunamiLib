const std = @import("std");
const contract = @import("../contract.zig");

/// Systematic Cauchy Reed-Solomon ERASURE code over GF(2^8) (poly 0x11D, α=2):
/// any `m` losses of `k+m` segments at KNOWN positions reconstruct bit-exact.
/// Erasure-only — it locates nothing; the caller marks which slots are missing.
pub const GF_POLY: u16 = 0x11D; // x^8 + x^4 + x^3 + x^2 + 1
pub const MAX_SEGS: usize = 256; // k + m ≤ 256 (distinct Cauchy nodes in GF(256))

// EXP doubled to 512 so mul can add logs without a %255 branch; LOG[0] is
// never read (the 0-guards below cover it).
const Tables = struct { exp: [512]u8, log: [256]u8 };
const TBL: Tables = blk: {
    @setEvalBranchQuota(20_000);
    var exp: [512]u8 = undefined;
    var log: [256]u8 = [_]u8{0} ** 256;
    var x: u16 = 1;
    var i: usize = 0;
    while (i < 255) : (i += 1) {
        exp[i] = @intCast(x);
        log[@intCast(x)] = @intCast(i);
        x <<= 1;
        if (x & 0x100 != 0) x ^= GF_POLY;
    }
    var j: usize = 255;
    while (j < 512) : (j += 1) exp[j] = exp[j - 255];
    break :blk .{ .exp = exp, .log = log };
};

pub const EXP: [512]u8 = TBL.exp;
pub const LOG: [256]u8 = TBL.log;

pub inline fn mul(a: u8, b: u8) u8 {
    if (a == 0 or b == 0) return 0;
    return EXP[@as(usize, LOG[a]) + @as(usize, LOG[b])];
}

/// Full product table: hot loops fetch `MUL[coeff]` once per row, then one lookup per byte.
pub const MUL: [256][256]u8 = blk: {
    @setEvalBranchQuota(1 << 20);
    var t: [256][256]u8 = undefined;
    for (0..256) |a| for (0..256) |b| {
        t[a][b] = mul(a, b);
    };
    break :blk t;
};

pub inline fn inv(a: u8) u8 {
    std.debug.assert(a != 0);
    return EXP[255 - @as(usize, LOG[a])];
}

pub inline fn div(a: u8, b: u8) u8 {
    std.debug.assert(b != 0);
    if (a == 0) return 0;
    return EXP[@as(usize, LOG[a]) + 255 - @as(usize, LOG[b])];
}

pub const Error = error{
    TooManySegments, // k + m > MAX_SEGS
    TooManyErasures, // more than m segments erased — unrecoverable
    Singular, // should never happen for a Cauchy survivor set
};

/// Systematic Cauchy RS coder for one (k, m). Owns the m×k generator matrix.
pub const Coder = struct {
    k: usize,
    m: usize,
    /// G[p][d] = g[p*k + d] = 1 / (x_p ⊕ y_d), x_p = k+p, y_d = d.
    g: []u8,

    pub fn init(alloc: std.mem.Allocator, k: usize, m: usize) !Coder {
        if (k + m > MAX_SEGS) return Error.TooManySegments;
        const g = try alloc.alloc(u8, k * m);
        for (0..m) |p| {
            const xp: u8 = @intCast(k + p);
            for (0..k) |d| {
                const yd: u8 = @intCast(d);
                g[p * k + d] = inv(xp ^ yd); // xp ^ yd != 0 (disjoint node sets)
            }
        }
        return .{ .k = k, .m = m, .g = g };
    }

    pub fn deinit(self: *Coder, alloc: std.mem.Allocator) void {
        alloc.free(self.g);
        self.* = undefined;
    }

    /// parity[p][i] = Σ_d G[p][d]·data[d][i]. `data` has k slices, `parity`
    /// has m slices, all of length `stripe_len`.
    pub fn encode(self: Coder, data: []const []const u8, parity: []const []u8, stripe_len: usize) void {
        contract.require(data.len == self.k and parity.len == self.m, "rs.encode: slice counts != k, m");
        for (data) |d| contract.require(d.len >= stripe_len, "rs.encode: data slice < stripe_len");
        for (parity) |p| contract.require(p.len >= stripe_len, "rs.encode: parity slice < stripe_len");
        for (0..self.m) |p| {
            const row = self.g[p * self.k ..][0..self.k];
            const out = parity[p];
            @memset(out[0..stripe_len], 0);
            for (0..self.k) |d| {
                const coeff = row[d];
                if (coeff == 0) continue;
                const mrow = &MUL[coeff];
                for (out[0..stripe_len], data[d][0..stripe_len]) |*o, s| o.* ^= mrow[s];
            }
        }
    }

    /// Reconstruct erased DATA slots in place. `segments` is k+m slices
    /// (0..k data, k..k+m parity), each `stripe_len` bytes; `present[j]`
    /// marks slot j as a trustworthy known row. Erased data slices may hold
    /// garbage on entry and are overwritten; erased parity slices are left
    /// untouched (never needed downstream). Returns TooManyErasures if more
    /// than m slots are missing.
    pub fn reconstruct(
        self: Coder,
        alloc: std.mem.Allocator,
        segments: []const []u8,
        present: []const bool,
        stripe_len: usize,
    ) !void {
        const k = self.k;
        const n = k + self.m;
        contract.require(segments.len == n and present.len == n, "rs.reconstruct: slice counts != k+m");
        for (segments) |s| contract.require(s.len >= stripe_len, "rs.reconstruct: segment < stripe_len");

        var erased_data = false;
        var n_present: usize = 0;
        for (0..n) |j| {
            if (present[j]) n_present += 1 else if (j < k) erased_data = true;
        }
        if (n_present < k) return Error.TooManyErasures;
        if (!erased_data) return; // only parity lost — nothing to do

        // Choose k surviving slots (any k; survivors ≥ k). Record each row's
        // source: a surviving data slot d contributes unit row e_d; a
        // surviving parity slot k+p contributes G[p][·].
        const surv = try alloc.alloc(usize, k);
        defer alloc.free(surv);
        {
            var r: usize = 0;
            var j: usize = 0;
            while (r < k) : (j += 1) {
                if (present[j]) {
                    surv[r] = j;
                    r += 1;
                }
            }
        }

        // A (k×k): row r = Gsys[surv[r]].
        const a = try alloc.alloc(u8, k * k);
        defer alloc.free(a);
        @memset(a, 0);
        for (0..k) |r| {
            const s = surv[r];
            if (s < k) {
                a[r * k + s] = 1; // unit row
            } else {
                const p = s - k;
                @memcpy(a[r * k ..][0..k], self.g[p * k ..][0..k]);
            }
        }

        // Ainv = A^{-1} over GF(256), Gauss-Jordan with partial pivoting.
        const ainv = try alloc.alloc(u8, k * k);
        defer alloc.free(ainv);
        try invert(alloc, a, ainv, k);

        // data[d] = Σ_r Ainv[d][r]·survivor[r], streamed row by row so each
        // pass is a linear table-lookup XOR. Survivors never alias erased slots.
        for (0..k) |d| {
            if (present[d]) continue;
            const dst = segments[d][0..stripe_len];
            @memset(dst, 0);
            for (ainv[d * k ..][0..k], 0..) |coeff, r| {
                if (coeff == 0) continue;
                const mrow = &MUL[coeff];
                for (dst, segments[surv[r]][0..stripe_len]) |*o, s| o.* ^= mrow[s];
            }
        }
    }
};

/// Invert k×k matrix `a` (row-major, consumed) into `out` over GF(256).
fn invert(alloc: std.mem.Allocator, a: []u8, out: []u8, k: usize) !void {
    // Augment [a | I] in a scratch k×2k buffer.
    const w = 2 * k;
    const aug = try alloc.alloc(u8, k * w);
    defer alloc.free(aug);
    @memset(aug, 0);
    for (0..k) |r| {
        @memcpy(aug[r * w ..][0..k], a[r * k ..][0..k]);
        aug[r * w + k + r] = 1;
    }
    for (0..k) |col| {
        // Partial pivot: Cauchy guarantees a nonzero pivot EXISTS (MDS), but
        // the mixed identity+Cauchy row order does not guarantee nonzero
        // leading minors, so a swap-in is required for correctness.
        if (aug[col * w + col] == 0) {
            var found = false;
            for (col + 1..k) |r| {
                if (aug[r * w + col] != 0) {
                    for (0..w) |j| {
                        const t = aug[col * w + j];
                        aug[col * w + j] = aug[r * w + j];
                        aug[r * w + j] = t;
                    }
                    found = true;
                    break;
                }
            }
            if (!found) return Error.Singular;
        }
        const pinv = inv(aug[col * w + col]);
        for (0..w) |j| aug[col * w + j] = mul(aug[col * w + j], pinv);
        for (0..k) |r| {
            if (r == col) continue;
            const f = aug[r * w + col];
            if (f == 0) continue;
            for (0..w) |j| aug[r * w + j] ^= mul(f, aug[col * w + j]);
        }
    }
    for (0..k) |r| @memcpy(out[r * k ..][0..k], aug[r * w + k ..][0..k]);
}

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;

test "field axioms" {
    try testing.expectEqual(@as(u8, 1), EXP[255]);
    try testing.expectEqual(@as(u8, 1), EXP[0]);
    var a: u16 = 1;
    while (a < 256) : (a += 1) {
        const av: u8 = @intCast(a);
        try testing.expectEqual(av, EXP[LOG[av]]); // EXP∘LOG identity
        try testing.expectEqual(@as(u8, 1), mul(av, inv(av))); // a·a⁻¹ = 1
        try testing.expectEqual(@as(u8, 0), av ^ av); // a + a = 0
        try testing.expectEqual(av, mul(av, 1)); // a·1 = a
        try testing.expectEqual(@as(u8, 0), mul(av, 0)); // a·0 = 0
        var b: u16 = 1;
        while (b < 256) : (b += 1) {
            const bv: u8 = @intCast(b);
            try testing.expectEqual(av, div(mul(av, bv), bv)); // (a·b)/b = a
        }
    }
}

test Coder {
    const alloc = testing.allocator;
    const k = 4;
    const m = 2;
    const L = 7; // stripe length
    var coder = try Coder.init(alloc, k, m);
    defer coder.deinit(alloc);

    // deterministic data
    var store: [k + m][L]u8 = undefined;
    var prng = std.Random.DefaultPrng.init(0x6A5);
    const rand = prng.random();
    for (0..k) |d| for (0..L) |i| {
        store[d][i] = rand.int(u8);
    };
    var data: [k][]const u8 = undefined;
    for (0..k) |d| data[d] = &store[d];
    var parity: [m][]u8 = undefined;
    for (0..m) |p| parity[p] = &store[k + p];
    coder.encode(&data, &parity, L);

    // Save a pristine copy to compare against after reconstruction.
    var orig: [k + m][L]u8 = store;

    // Enumerate every erasure pattern of size 0, 1, 2 over k+m slots.
    const n = k + m;
    for (0..n + 1) |e1| {
        for (e1..n + 1) |e2_| {
            // patterns: {} (e1==n flag), {e1}, {e1,e2}
            var present = [_]bool{true} ** n;
            var n_er: usize = 0;
            if (e1 < n) {
                present[e1] = false;
                n_er += 1;
            }
            if (e2_ < n and e2_ != e1) {
                present[e2_] = false;
                n_er += 1;
            }
            if (n_er > m) continue;
            // restore garbage into erased slots to prove reconstruction
            var work: [k + m][L]u8 = orig;
            for (0..n) |j| if (!present[j]) {
                @memset(&work[j], 0xAA);
            };
            var segs: [k + m][]u8 = undefined;
            for (0..n) |j| segs[j] = &work[j];
            try coder.reconstruct(alloc, &segs, &present, L);
            // every data slot must match the original
            for (0..k) |d| try testing.expectEqualSlices(u8, &orig[d], &work[d]);
        }
    }
}

test "MDS pin: every size-2 erasure of k=4,m=2 is invertible & exact" {
    // Covered by the round-trip test's exhaustive enumeration above; this
    // pin states the property explicitly for the record.
    const alloc = testing.allocator;
    var coder = try Coder.init(alloc, 4, 2);
    defer coder.deinit(alloc);
    const n = 6;
    for (0..n) |i| {
        for (i + 1..n) |j| {
            var present = [_]bool{true} ** n;
            present[i] = false;
            present[j] = false;
            // build the survivor matrix and confirm it inverts
            var surv: [4]usize = undefined;
            var r: usize = 0;
            var s: usize = 0;
            while (r < 4) : (s += 1) {
                if (present[s]) {
                    surv[r] = s;
                    r += 1;
                }
            }
            var a = [_]u8{0} ** 16;
            for (0..4) |rr| {
                const sl = surv[rr];
                if (sl < 4) a[rr * 4 + sl] = 1 else @memcpy(a[rr * 4 ..][0..4], coder.g[(sl - 4) * 4 ..][0..4]);
            }
            var ainv = [_]u8{0} ** 16;
            try invert(alloc, &a, &ainv, 4);
        }
    }
}

test "too many erasures fails gracefully" {
    const alloc = testing.allocator;
    const k = 4;
    const m = 2;
    const L = 4;
    var coder = try Coder.init(alloc, k, m);
    defer coder.deinit(alloc);
    var work: [k + m][L]u8 = undefined;
    @memset(std.mem.asBytes(&work), 0);
    var segs: [k + m][]u8 = undefined;
    for (0..k + m) |j| segs[j] = &work[j];
    // erase 3 (one parity + two data) > m=2
    var present = [_]bool{true} ** (k + m);
    present[0] = false;
    present[1] = false;
    present[k] = false;
    try testing.expectError(Error.TooManyErasures, coder.reconstruct(alloc, &segs, &present, L));
}
