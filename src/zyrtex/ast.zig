const std = @import("std");
const Array = std.ArrayList;
const Allocator = std.mem.Allocator;
const tst = std.testing;
const math = std.math;

pub const GenericNode = struct { kind: NodeKind, renderInfo: ?RenderInfo = null, position: ?Position = null };
pub const NodeKind = enum {
    root,
    string,
    whitespace,
    parbreak,
    comment,
    macro,
    environment,
    verbatimEnvironment,
    inlineMath,
    displayMath,
    group,
    verb,
};
pub const RenderInfo = struct {
    inMathMode: bool = false,
    alignContent: bool = false,
    pgfkeysArgs: bool = false,
    tikzEnvironment: bool = false,
    tikzPathComment: bool = false,
    namedArguments: ?[][]const u8 = null,
    inParMode: bool = false,
    breakAfter: bool = false,
    breakAround: bool = false,
    breakBefore: bool = false,
    hangingIndent: bool = false,
    sysdelims: ?[]Node = null,
};

pub const Position = struct {
    offset: usize = 0,
    line: usize = 0,
    column: usize = 0,
};
pub const Macro = struct {
    content: []const u8,
    escapeToken: ?[]const u8 = null,
    args: ?[]Argument = null,
};
pub const Comment = struct {
    content: []const u8,
    sameLine: bool = false,
    suffixParbreak: bool = false,
    leadingWhitespace: bool = false,
};
pub const Environment = struct {
    isMath: bool,
    env: []const u8,
    args: ?[]Argument = null,
};
pub const VerbatimEnvironment = struct {
    env: []const u8,
    content: []const u8,
};

pub const DisplayMath = struct {};
pub const InlineMath = struct {};
pub const Group = struct {};
pub const ParBreak = struct {};
pub const Root = struct {};
pub const Whitespace = struct {};
pub const String = struct {
    content: []const u8,
};
pub const Node = union(NodeKind) {
    verb: Verb,
    group: Group,
    displayMath: DisplayMath,
    inlineMath: InlineMath,
    environment: Environment,
    verbatimEnvironment: VerbatimEnvironment,
    macro: Macro,
    comment: Comment,
    parbreak: ParBreak,
    whitespace: Whitespace,
    string: String,
    root: Root,
};

pub const Verb = struct {
    env: []const u8,
    escape: []const u8,
    content: []const u8,
};

pub const Argument = struct {
    openMark: []const u8,
    closeMark: []const u8,
};
pub const EnvInfo = struct {
    renderInfo: ?RenderInfo = null,
    processContent: ?*const fn ([]Node) anyerror![]Node = null,
    signature: ?[]const u8 = null,
};
pub const MacroInfo = struct {
    renderInfo: ?RenderInfo = null,
    signature: ?[]const u8 = null,
    escapeToken: ?[]const u8 = null,
    argumentParser: ?ArgumentParser = null,
};
pub const ArgumentParserResult = struct {
    arguments: []Argument,
    nodesRemoved: usize = 0,
};
pub const ArgumentParser = *const fn (node: []Node, startPos: usize) anyerror!ArgumentParserResult;
