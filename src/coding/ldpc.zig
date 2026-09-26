const std = @import("std");
pub const tables = @import("ldpc_tables.zig");

/// Normalized min-sum scale: raw min-sum overestimates check magnitudes by
/// ~0.3–0.5 dB at row degree ~20; α≈0.8 lands within ~0.1 dB of sum-product.
pub const minsum_alpha: f32 = 0.8;
pub const default_iters: u32 = 50;

/// Quasi-cyclic LDPC with the 802.11n parity structure (weight-3 h-column +
/// shift-0 dual diagonal), so encoding is a linear recursion rather than
/// Gaussian elimination. Edge tables are expanded at comptime.
///
/// LLR convention: positive = bit more likely 1. The check-node sign flip
/// below is coupled to it.
pub fn Code(comptime Z: usize, comptime base: anytype) type {
    const mb = base.len;
    const nb = base[0].len;
    comptime std.debug.assert(mb >= 3 and nb > mb);
    return struct {
        pub const n: usize = nb * Z;
        pub const m: usize = mb * Z;
        pub const k: usize = n - m;
        const kb: usize = nb - mb;
        const Var = std.math.IntFittingRange(0, n - 1);

        pub const num_blocks: usize = blk: {
            var cnt: usize = 0;
            for (base) |row| for (row) |v| {
                cnt += @intFromBool(v >= 0);
            };
            break :blk cnt;
        };
        pub const num_edges: usize = num_blocks * Z;

        pub const row_deg: [mb]usize = blk: {
            var deg: [mb]usize = undefined;
            for (base, 0..) |row, i| {
                var cnt: usize = 0;
                for (row) |v| cnt += @intFromBool(v >= 0);
                deg[i] = cnt;
            }
            break :blk deg;
        };
        const max_deg = std.mem.max(usize, &row_deg);

        /// Ordered by block row, expanded row z, then non-null block column:
        /// row (i, z) touches variable j·Z + ((z + s_ij) mod Z).
        const edge_var: [num_edges]Var = blk: {
            @setEvalBranchQuota(num_edges * 64);
            var ev: [num_edges]Var = undefined;
            var e: usize = 0;
            for (base) |row| for (0..Z) |z| for (row, 0..) |s, j| {
                if (s >= 0) {
                    ev[e] = j * Z + ((z + s) % Z);
                    e += 1;
                }
            };
            break :blk ev;
        };

        /// q0 = Σλ_i (h-column top/bottom shifts cancel, the interior shift-0
        /// survives); q_{i+1} = q_i + λ_i + rot(q0, h_i).
        pub fn encode(info: *const [k]u1, out: *[n]u1) void {
            @memcpy(out[0..k], info);

            var lambda: [mb][Z]u1 = @splat(@splat(0));
            inline for (0..mb) |i| inline for (0..kb) |j| {
                const s = base[i][j];
                if (s >= 0) {
                    const blk_in = info[j * Z ..][0..Z];
                    for (0..Z) |r| lambda[i][r] ^= blk_in[(r + s) % Z];
                }
            };

            var q0: [Z]u1 = @splat(0);
            for (lambda) |l| for (&q0, l) |*q, v| {
                q.* ^= v;
            };
            out[k..][0..Z].* = q0;

            const h_top: usize = base[0][kb];
            var q_prev: [Z]u1 = undefined;
            for (0..Z) |r| q_prev[r] = lambda[0][r] ^ q0[(r + h_top) % Z];
            out[k + Z ..][0..Z].* = q_prev;

            inline for (1..mb - 1) |i| {
                const h = base[i][kb];
                for (0..Z) |r| {
                    var v = q_prev[r] ^ lambda[i][r];
                    if (h >= 0) v ^= q0[(r + h) % Z];
                    q_prev[r] = v;
                }
                out[k + (i + 1) * Z ..][0..Z].* = q_prev;
            }

            // Last block row is redundant — it must balance.
            if (std.debug.runtime_safety) {
                const h_bot: usize = base[mb - 1][kb];
                for (0..Z) |r| std.debug.assert(lambda[mb - 1][r] ^ q0[(r + h_bot) % Z] ^ q_prev[r] == 0);
            }
        }

        /// Walks the expanded edge tables — independent of the encoder's
        /// recursion, so it catches expansion bugs.
        pub fn checkParity(cw: *const [n]u1) bool {
            var e: usize = 0;
            inline for (row_deg) |deg| for (0..Z) |_| {
                var x: u1 = 0;
                inline for (0..deg) |d| x ^= cw[edge_var[e + d]];
                if (x != 0) return false;
                e += deg;
            };
            return true;
        }

        /// Flooding normalized min-sum BP. Iteration 0 checks raw hard
        /// decisions. Returns true when all checks pass; `out` always holds
        /// the current decisions.
        pub fn decode(llr: *const [n]f32, out: *[n]u1, max_iters: u32) bool {
            var r_msg: [num_edges]f32 = @splat(0);
            var total: [n]f32 = undefined;

            var iter: u32 = 0;
            while (true) : (iter += 1) {
                total = llr.*;
                for (r_msg, edge_var) |rv, v| total[v] += rv;
                for (out, total) |*b, t| b.* = @intFromBool(t > 0);
                if (checkParity(out)) return true;
                if (iter == max_iters) return false;

                // Extrinsic sign under positive=1 is (−1)^deg · Π sign(others),
                // so the flip depends on each block row's degree parity.
                var e: usize = 0;
                inline for (row_deg) |deg| {
                    const deg_sign: f32 = if (deg % 2 == 1) -minsum_alpha else minsum_alpha;
                    for (0..Z) |_| {
                        var q: [max_deg]f32 = undefined;
                        var min1: f32 = std.math.floatMax(f32);
                        var min2: f32 = std.math.floatMax(f32);
                        var neg: bool = false;
                        inline for (0..deg) |d| {
                            const qv = total[edge_var[e + d]] - r_msg[e + d];
                            q[d] = qv;
                            const a = @abs(qv);
                            if (a < min1) {
                                min2 = min1;
                                min1 = a;
                            } else if (a < min2) min2 = a;
                            neg = neg != (qv < 0);
                        }
                        inline for (0..deg) |d| {
                            const mag = if (@abs(q[d]) == min1) min2 else min1;
                            const flip = neg != (q[d] < 0);
                            r_msg[e + d] = if (flip) -deg_sign * mag else deg_sign * mag;
                        }
                        e += deg;
                    }
                }
            }
        }
    };
}

pub const Ldpc648R12 = Code(27, tables.BASE_648_R12);
pub const Ldpc1944R12 = Code(81, tables.BASE_1944_R12);
pub const Ldpc1944R23 = Code(81, tables.BASE_1944_R23);
pub const Ldpc1944R34 = Code(81, tables.BASE_1944_R34);
pub const Ldpc1944R56 = Code(81, tables.BASE_1944_R56);

// ── Tests ────────────────────────────────────────────────────────────────

const testing = std.testing;
const golden = @import("ldpc_golden.zig");

fn checkStructure(comptime base: anytype, comptime Z: usize) !void {
    const mb = base.len;
    const kb = 24 - mb;
    for (base) |row| for (row) |s| try testing.expect(s >= -1 and s < Z);
    try testing.expectEqual(base[0][kb], base[mb - 1][kb]);
    var interior: usize = 0;
    for (1..mb - 1) |i| if (base[i][kb] >= 0) {
        try testing.expectEqual(@as(i16, 0), base[i][kb]);
        interior += 1;
    };
    try testing.expectEqual(@as(usize, 1), interior);
    for (0..mb - 1) |t| for (0..mb) |i| {
        const want: i16 = if (i == t or i == t + 1) 0 else -1;
        try testing.expectEqual(want, base[i][kb + 1 + t]);
    };
}

test "base matrix structure" {
    try checkStructure(tables.BASE_1944_R12, 81);
    try checkStructure(tables.BASE_1944_R23, 81);
    try checkStructure(tables.BASE_1944_R34, 81);
    try checkStructure(tables.BASE_1944_R56, 81);
    try checkStructure(tables.BASE_648_R12, 27);
    try testing.expectEqualSlices(usize, &.{ 20, 20, 20, 19 }, &Ldpc1944R56.row_deg);
    const total = Ldpc648R12.num_blocks + Ldpc1944R12.num_blocks +
        Ldpc1944R23.num_blocks + Ldpc1944R34.num_blocks + Ldpc1944R56.num_blocks;
    try testing.expectEqual(@as(usize, 426), total);
}

fn goldenCheck(comptime C: type, infos: anytype, cws: anytype) !void {
    for (infos, cws) |info, want| {
        var got: [C.n]u1 = undefined;
        C.encode(&info, &got);
        try testing.expectEqualSlices(u1, &want, &got);
        try testing.expect(C.checkParity(&got));
    }
}

test "golden encode vectors" {
    try goldenCheck(Ldpc1944R12, golden.GOLDEN_1944_R12_INFO, golden.GOLDEN_1944_R12_CW);
    try goldenCheck(Ldpc1944R23, golden.GOLDEN_1944_R23_INFO, golden.GOLDEN_1944_R23_CW);
    try goldenCheck(Ldpc1944R34, golden.GOLDEN_1944_R34_INFO, golden.GOLDEN_1944_R34_CW);
    try goldenCheck(Ldpc1944R56, golden.GOLDEN_1944_R56_INFO, golden.GOLDEN_1944_R56_CW);
    try goldenCheck(Ldpc648R12, golden.GOLDEN_648_R12_INFO, golden.GOLDEN_648_R12_CW);
}

test "decode: clean, sign flips, erasure burst (r5/6)" {
    const C = Ldpc1944R56;
    const cw = golden.GOLDEN_1944_R56_CW[2];
    var llr: [C.n]f32 = undefined;
    for (&llr, cw) |*l, b| l.* = if (b == 1) 4.0 else -4.0;
    var out: [C.n]u1 = undefined;

    try testing.expect(C.decode(&llr, &out, 0));
    try testing.expectEqualSlices(u1, &cw, &out);

    var flipped = llr;
    var i: usize = 13;
    for (0..15) |_| {
        flipped[i] = -flipped[i];
        i += 127;
    }
    try testing.expect(C.decode(&flipped, &out, default_iters));
    try testing.expectEqualSlices(u1, &cw, &out);

    var erased = llr;
    @memset(erased[700..740], 0);
    try testing.expect(C.decode(&erased, &out, default_iters));
    try testing.expectEqualSlices(u1, &cw, &out);
}

test "waterfall: pre-FEC ~1e-2 decodes near-clean (r5/6)" {
    const C = Ldpc1944R56;
    var prng: std.Random.DefaultPrng = .init(0xFEC);
    const r = prng.random();
    const sigma: f32 = 0.4299;
    var info: [C.k]u1 = undefined;
    var cw: [C.n]u1 = undefined;
    var llr: [C.n]f32 = undefined;
    var out: [C.n]u1 = undefined;
    var fails: usize = 0;
    for (0..100) |_| {
        for (&info) |*b| b.* = r.int(u1);
        C.encode(&info, &cw);
        for (&llr, cw) |*l, b| {
            const y = @as(f32, if (b == 1) 1.0 else -1.0) + sigma * r.floatNorm(f32);
            l.* = 2.0 * y / (sigma * sigma);
        }
        const ok = C.decode(&llr, &out, default_iters);
        if (!ok or !std.mem.eql(u1, &out, &cw)) fails += 1;
    }
    try testing.expect(fails <= 2);
}
