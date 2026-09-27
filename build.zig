const std = @import("std");

/// `tsunami` is std-only. `tsunami_rl` exists only with `-Drl=true`, so the
/// lazy raylib dependency is never fetched by consumers that don't ask for it.
pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const rl = b.option(bool, "rl", "Expose the tsunami_rl raylib module") orelse false;
    const test_filters = b.option(
        []const []const u8,
        "test-filter",
        "Skip tests that do not match any of the specified filters",
    ) orelse &.{};

    // ── Modules ──────────────────────────────────────────────────────────
    const core = b.addModule("tsunami", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    if (rl) _ = rlModule(b, core, target, optimize, true);

    // ── Tests ────────────────────────────────────────────────────────────
    const test_step = b.step("test", "Run unit tests");
    addSuite(b, test_step, target, optimize, test_filters, rl);

    // ── Examples ─────────────────────────────────────────────────────────
    // `examples` builds and runs every one, so they can't rot; `run-<name>` runs one.
    const examples_step = b.step("examples", "Build and run every example");
    for (examples) |name| {
        const m = b.createModule(.{
            .root_source_file = b.path(b.fmt("examples/{s}.zig", .{name})),
            .target = target,
            .optimize = optimize,
        });
        m.addImport("tsunami", core);
        const exe = b.addExecutable(.{ .name = name, .root_module = m });
        const run = b.addRunArtifact(exe);
        if (b.args) |a| run.addArgs(a);
        b.step(b.fmt("run-{s}", .{name}), b.fmt("Run examples/{s}.zig", .{name})).dependOn(&run.step);
        examples_step.dependOn(&run.step);
    }

    // ── Docs ─────────────────────────────────────────────────────────────
    const docs_obj = b.addObject(.{ .name = "tsunami", .root_module = core });
    const docs = b.addInstallDirectory(.{
        .source_dir = docs_obj.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    });
    b.step("docs", "Emit autodoc HTML to zig-out/docs").dependOn(&docs.step);

    const serve_exe = b.addExecutable(.{ .name = "docs_serve", .root_module = b.createModule(.{
        .root_source_file = b.path("tools/docs_serve.zig"),
        .target = target,
        .optimize = optimize,
    }) });
    const serve = b.addRunArtifact(serve_exe);
    serve.addArg(b.getInstallPath(.prefix, "docs"));
    if (b.args) |a| serve.addArgs(a);
    serve.step.dependOn(&docs.step);
    b.step("docs-serve", "Build docs and serve on 127.0.0.1:8080 (port override: -- <port>)").dependOn(&serve.step);

    // ── CI gate ──────────────────────────────────────────────────────────
    const ci_step = b.step("ci", "Format check + tests in Debug/ReleaseSafe/ReleaseFast");
    const fmt = b.addFmt(.{ .paths = &.{ "src", "examples", "tools", "build.zig", "build.zig.zon" }, .check = true });
    ci_step.dependOn(&fmt.step);
    ci_step.dependOn(examples_step);
    for ([_]std.builtin.OptimizeMode{ .Debug, .ReleaseSafe, .ReleaseFast }) |mode| {
        addSuite(b, ci_step, target, mode, test_filters, rl);
    }
}

fn rlModule(
    b: *std.Build,
    core: *std.Build.Module,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    exported: bool,
) ?*std.Build.Module {
    const dep = b.lazyDependency("raylib_zig", .{ .target = target, .optimize = optimize }) orelse return null;
    const opts: std.Build.Module.CreateOptions = .{
        .root_source_file = b.path("src/rl/root.zig"),
        .target = target,
        .optimize = optimize,
    };
    const m = if (exported) b.addModule("tsunami_rl", opts) else b.createModule(opts);
    m.addImport("tsunami", core);
    m.addImport("raylib", dep.module("raylib"));
    m.linkLibrary(dep.artifact("raylib"));
    return m;
}

fn addSuite(
    b: *std.Build,
    step: *std.Build.Step,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    filters: []const []const u8,
    rl: bool,
) void {
    const core = b.createModule(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    const t = b.addTest(.{ .name = "tsunami", .root_module = core, .filters = filters });
    step.dependOn(&b.addRunArtifact(t).step);

    if (!rl) return;
    const m = rlModule(b, core, target, optimize, false) orelse return;
    const rt = b.addTest(.{ .name = "tsunami_rl", .root_module = m, .filters = filters });
    step.dependOn(&b.addRunArtifact(rt).step);
}

const examples = [_][]const u8{
    "spectrum",
    "modem",
    "erasure",
    "geometry",
    "cli",
    "calculator",
    "wav",
    "cpu",
    "dataflow",
    "animation",
    "collections",
    "terminal",
};
