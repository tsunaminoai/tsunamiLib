//! Always-on API-entry precondition checks. `std.debug.assert` vanishes in
//! ReleaseFast, turning a caller bug (mismatched slice lengths) into silent
//! out-of-bounds; `require` survives every release mode.

const builtin = @import("builtin");

/// Once per API entry, never per element — keep `std.debug.assert` there.
pub inline fn require(ok: bool, comptime what: []const u8) void {
    if (!ok) {
        // 0.16's default panic handler doesn't compile for emscripten.
        if (comptime builtin.target.os.tag == .emscripten) @trap();
        @panic("tsunami: precondition violated: " ++ what);
    }
}
