# tsunamiLib

Zig 0.16.0. The core module `tsunami` depends only on std. `tsunami_rl` (raylib) is opt-in with `-Drl=true`, and raylib is a lazy dependency, so it is never fetched unless you ask for it.

```zig
// build.zig.zon: .tsunamiLib = .{ .url = "...", .hash = "..." }
const ts = b.dependency("tsunamiLib", .{ .target = target, .optimize = optimize });
mod.addImport("tsunami", ts.module("tsunami"));
// raylib helpers: b.dependency("tsunamiLib", .{ ..., .rl = true }).module("tsunami_rl")
```

`zig build test` · `zig build ci` (fmt check, tests in Debug/ReleaseSafe/ReleaseFast, all examples) · `tools/zt src/<file>.zig` (runs one file's tests)

**Docs:** `zig build docs-serve`, then open http://127.0.0.1:8080/. `zig build docs` only writes `zig-out/docs`, which has to be served over HTTP. Each major declaration's page includes a doctest example.

**Examples** (`zig build run-<name>`, `zig build examples` runs them all). Each one checks its own result and exits non-zero on failure.

| Example | Shows |
|---|---|
| `spectrum` | comptime FIR taps, SIMD `Fir`, biquad, windowed `Fft` |
| `modem` | 64-QAM → AWGN → LLRs → LDPC(1944, r5/6) decode |
| `erasure` | Reed-Solomon reconstruct after losing 2 of 6 segments; interleaver burst spreading |
| `geometry` | vec/Mat4 projection, ray–AABB, bezier distance, great-circle, grid, GMST |
| `wav` | WAV round trip through `std.Io`; save blob with corruption rejected |
| `cli` | struct-driven `args.parse`, comptime usage, diagnostics |
| `calculator` | comptime-spec lexer feeding the Pratt parser |
| `cpu` | comptime opcode `Isa` over a memory-mapped `Bus` (sums 1..10) |
| `dataflow` | node `Graph` evaluated in topological order, cycle rejection |
| `animation` | tweens, typewriter, countdown, scene stack, stopwatch |
| `collections` | seeded shuffle, event log, interner, typed ids |
| `terminal` | comptime ANSI, braille sparkline, terminal size |

Conventions:
- Sizes, tables and dispatch are comptime.
- Callers pass scratch buffers, and hot paths don't allocate.
- Where allocation is needed, the allocator is the first parameter.
- `contract.require` checks caller contracts and stays on in every release mode.
- Only `io/*` and `app/*` touch `std.Io`.

| Module | Provides | Origin |
|---|---|---|
| `math.vec` | `@Vector` ops, swizzle, `Mat(n,T)`, translate/rotate/perspective/lookAt | new |
| `math.grid` | `Grid(T)`, `Grid3`, `Pos(I)`, neighbour iterators | zephyris, aoc2024 |
| `math.geo` | great-circle, bearing, ENU, 4/3-earth radar beam | radar, zephyris |
| `math.geom` | ray/AABB, point–segment/bezier distance, `Rect` | radar, node-editor |
| `math.astro` | GMST, hour angle, spherical↔cartesian | tsunamiLib |
| `dsp.fft` | `Fft(T,n)` radix-2 with comptime twiddles/bitrev | tinyradio, soundtracer |
| `dsp.window` / `fir` | comptime windows, lowpass/RRC taps, SIMD `Fir(T,taps)` | tapes, tinyradio |
| `dsp.biquad` | RBJ `Biquad(T)`, `Cascade`, Butterworth | tapes (legacy, rewritten) |
| `dsp.loop` | PI `LoopFilter(T)`, `Pll(T)` | tapes |
| `dsp.interp` / `resample` | polyphase sinc, Farrow, rational resampler `Plan` | tapes |
| `dsp.sample` / `channel` | PCM conversion; AWGN, wow/flutter, dropout | tapes |
| `coding.ldpc` | QC-LDPC `Code(Z, base)`, 802.11n tables, min-sum decode | tapes |
| `coding.rs` | Cauchy RS erasure `Coder` over GF(256) | tapes |
| `coding.constellation` | Gray `Pam`/`Qam` with O(1) slicing and exact max-log LLRs | tapes (rewritten) |
| `coding.interleave` | block (de)interleaver | tapes |
| `collections.*` | `Pile(T)`, `EventLog(E)`, typed `Id`/`IdSlab`, string `Interner` | cribbage, zero |
| `rand.rng` | replayable seeded `Rng` | cribbage |
| `text.args` | struct-driven CLI parser with comptime usage | the-wall |
| `text.lex` / `pratt` | comptime-spec zero-alloc lexer; depth-capped Pratt `Parser` | zero, monkey |
| `io.wav` / `blocks` / `blob` / `save` | RIFF over `std.Io`; extern-struct endian views; CRC'd save blobs; native/emscripten persistence | tapes, lastchoice, the-wall |
| `term.*` | TIOCGWINSZ size, comptime ANSI, braille | termui, cribbage |
| `app.*` | `HotReload` (shadow-copied dlopen), `SceneStack`, `Stopwatch`/`Countdown`, `Typewriter`, `Tween` | cribbage, trickdraft, the-wall |
| `emu.dispatch` / `bus` | comptime opcode `Isa` with threaded `run`; vtable memory `Bus` | nes |
| `graph.node` | tagged-union dataflow `Graph(V, kinds)` with static dispatch | node-editor |
| `tsunami_rl` | camera follow, theme, layout regions, draw helpers, headless screenshots | trickdraft, the-wall |

The raylib module needs X11/GL headers; `devbox shell` provides them.
