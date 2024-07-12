const std = @import("std");
const builtin = @import("builtin");
const print = std.debug.print;
const TermSize = @This();

width: usize,
height: usize,

pub fn init() !TermSize {
    return switch (builtin.os.tag) {
        .windows => {
            return undefined;
        },

        .macos, .linux, .freebsd => {
            return undefined;
        },
        else => {
            return undefined;
        },
    };
}

fn get_width_from_IO() !usize {
    const stdout = std.io.getStdOut();
    const s = try stdout.getEndPos();
    return s;
}

const unix_term = extern struct {
    ws_row: u16,
    ws_col: u16,
    ws_xpixel: u16,
    ws_ypixel: u16,

    const self = @This();
    pub fn init() !TermSize {
        var ws: self = undefined;
        const TIOCGWINSZ = 0x40087468;

        const fd = std.posix.STDOUT_FILENO;
        const rc = std.posix.system.ioctl(fd, TIOCGWINSZ, @intFromPtr(&ws));
        if (rc != 0) {
            std.debug.print("ioctl failed\n", .{});
            return error.IoCtlError;
        }
        return .{
            .width = ws.ws_col,
            .height = ws.ws_row,
        };
    }
};

test {
    const ts = try TermSize.init();
    print("{any}\n", .{ts});
}
