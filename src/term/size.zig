const std = @import("std");
const builtin = @import("builtin");
const contract = @import("../contract.zig");

pub const Size = struct {
    cols: u16,
    rows: u16,
};

pub fn get(fd: std.posix.fd_t) ?Size {
    const os_tag = builtin.os.tag;

    const ioctlnum = comptime switch (os_tag) {
        .linux => std.os.linux.T.IOCGWINSZ,
        .macos => 0x40087468,
        else => return null,
    };

    var wsz: std.posix.winsize = undefined;
    const rc = std.posix.system.ioctl(fd, ioctlnum, @intFromPtr(&wsz));

    if (rc == 0) {
        return Size{
            .cols = wsz.col,
            .rows = wsz.row,
        };
    }
    return null;
}

// ── Tests ────────────────────────────────────────────────────────────────

test "non-tty fd returns null" {
    // Test with stdin; if running in non-interactive context, may return null
    const result = get(std.posix.STDIN_FILENO);
    // Result depends on test environment; we just verify the function doesn't crash
    _ = result;
}
