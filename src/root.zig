//! tsunamiLib: std-only building blocks for Zig 0.16.0.
//!
//! Sizes, tables and dispatch are comptime; callers pass scratch buffers;
//! allocating APIs take the allocator first; only `io` and `app` touch
//! `std.Io`. Each major declaration's page shows a runnable doctest.
//! Full programs live in `examples/` (`zig build run-<name>`).
//! Raylib helpers are a separate module, `tsunami_rl` (`-Drl=true`).

const std = @import("std");

pub const contract = @import("contract.zig");

pub const math = struct {
    pub const vec = @import("math/vec.zig");
    pub const grid = @import("math/grid.zig");
    pub const geo = @import("math/geo.zig");
    pub const geom = @import("math/geom.zig");
    pub const astro = @import("math/astro.zig");
};

pub const dsp = struct {
    pub const fft = @import("dsp/fft.zig");
    pub const window = @import("dsp/window.zig");
    pub const fir = @import("dsp/fir.zig");
    pub const biquad = @import("dsp/biquad.zig");
    pub const loop = @import("dsp/loop.zig");
    pub const interp = @import("dsp/interp.zig");
    pub const resample = @import("dsp/resample.zig");
    pub const sample = @import("dsp/sample.zig");
    pub const channel = @import("dsp/channel.zig");
};

pub const coding = struct {
    pub const rs = @import("coding/rs.zig");
    pub const ldpc = @import("coding/ldpc.zig");
    pub const constellation = @import("coding/constellation.zig");
    pub const interleave = @import("coding/interleave.zig");
};

pub const collections = struct {
    pub const pile = @import("collections/pile.zig");
    pub const event_log = @import("collections/event_log.zig");
    pub const ids = @import("collections/ids.zig");
    pub const intern = @import("collections/intern.zig");
};

pub const rand = struct {
    pub const rng = @import("rand/rng.zig");
};

pub const text = struct {
    pub const args = @import("text/args.zig");
    pub const lex = @import("text/lex.zig");
    pub const pratt = @import("text/pratt.zig");
};

pub const io = struct {
    pub const wav = @import("io/wav.zig");
    pub const blocks = @import("io/blocks.zig");
    pub const blob = @import("io/blob.zig");
    pub const save = @import("io/save.zig");
};

pub const term = struct {
    pub const size = @import("term/size.zig");
    pub const ansi = @import("term/ansi.zig");
    pub const braille = @import("term/braille.zig");
};

pub const app = struct {
    pub const hotreload = @import("app/hotreload.zig");
    pub const scene = @import("app/scene.zig");
    pub const stopwatch = @import("app/stopwatch.zig");
    pub const typewriter = @import("app/typewriter.zig");
    pub const tween = @import("app/tween.zig");
};

pub const emu = struct {
    pub const dispatch = @import("emu/dispatch.zig");
    pub const bus = @import("emu/bus.zig");
};

pub const graph = struct {
    pub const node = @import("graph/node.zig");
};

test {
    inline for (comptime std.meta.declarations(@This())) |d| {
        const v = @field(@This(), d.name);
        if (@TypeOf(v) == type) std.testing.refAllDecls(v);
    }
}
