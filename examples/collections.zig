//! `collections.pile.Pile` shuffled deterministically via `rand.rng.Rng`,
//! `collections.intern.Interner` dedup, typed `collections.ids.Id`, and
//! `collections.event_log.EventLog` drain.
const std = @import("std");
const ts = @import("tsunami");
const collections = ts.collections;

// Two tags produce distinct types — mixing them up is a compile error, not a
// runtime bug: an ItemId can't be passed where a PlayerId is expected.
const ItemId = collections.ids.Id("item");
const PlayerId = collections.ids.Id("player");
comptime {
    std.debug.assert(ItemId != PlayerId);
}

const Event = union(enum) { drew: u8, discarded: u8 };

pub fn main(init: std.process.Init) !void {
    var buf: [1024]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(init.io, &buf);
    const w = &stdout.interface;

    var gpa_state: std.heap.DebugAllocator(.{}) = .init;
    defer std.debug.assert(gpa_state.deinit() == .ok);
    const gpa = gpa_state.allocator();

    var pile: collections.pile.Pile(u8) = .empty;
    defer pile.deinit(gpa);
    for (0..13) |i| try pile.append(gpa, @intCast(i));

    var rng = ts.rand.rng.Rng.init(0xC0FFEE);
    pile.shuffle(rng.random());
    var rng2 = ts.rand.rng.Rng.init(0xC0FFEE);
    var expect: collections.pile.Pile(u8) = .empty;
    defer expect.deinit(gpa);
    for (0..13) |i| try expect.append(gpa, @intCast(i));
    expect.shuffle(rng2.random());
    try w.print("shuffled: {any}\n", .{pile.items()});
    if (!std.mem.eql(u8, pile.items(), expect.items())) return error.ExampleFailed;

    var log: collections.event_log.EventLog(Event) = .empty;
    defer log.deinit(gpa);
    var i: usize = 0;
    while (pile.draw()) |card| : (i += 1) {
        try log.append(gpa, .{ .drew = card });
        if (i >= 2) break;
    }
    try w.print("drained {d} events: first={d}\n", .{ log.len(), log.drain()[0].drew });
    if (log.len() != 3) return error.ExampleFailed;

    var itn: collections.intern.Interner = .empty;
    defer itn.deinit(gpa);
    const a = try itn.intern(gpa, "sword");
    const b = try itn.intern(gpa, "shield");
    const c = try itn.intern(gpa, "sword"); // dedups to the same handle as `a`
    try w.print("intern: sword=={s} shield!=sword: {}\n", .{ itn.get(a), a != b });
    if (a != c or a == b) return error.ExampleFailed;

    const hero: PlayerId = .from(0);
    const sword: ItemId = .from(1);
    try w.print("ids: {f} {f}\n", .{ hero, sword });
    try w.flush();
}
