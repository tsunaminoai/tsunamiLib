const std = @import("std");
const assert = std.debug.assert;
const Array = std.ArrayList;

const NodeType = enum {
    Document,
    Section,
    Heading,
    Text,
    Bold,
    Italic,
    Code,
    ListItem,
    Paragraph,
    Link,
    Image,
    BlockQuote,
    HorizontalRule,
};

pub const Node = struct {
    ptr: *anyopaque,
    vtab: *const VTab,
    type: NodeType,

    const VTab = struct {
        renderFn: *const fn (*anyopaque, void) void,
        addChildFn: *const fn (*anyopaque, Node) void,
        deinitFn: *const fn (*anyopaque) void,
    };
    fn render(self: *Node, renderer: void) void {
        self.vtab.renderFn(renderer);
    }
    fn addChild(self: *Node, child: Node) void {
        self.vtab.addChildFn(child);
    }
    pub fn deinit(self: *Node) void {
        self.vtab.deinitFn(self);
    }

    pub fn init(obj: anytype) Node {
        const Ptr = @TypeOf(obj);
        const PtrInfo = @typeInfo(Ptr);
        assert(PtrInfo == .Pointer);
        assert(PtrInfo.Pointer.size == .One);
        assert(@typeInfo(PtrInfo.Pointer.child) == .Struct);

        const alignment = PtrInfo.Pointer.alignment;
        _ = alignment; // autofix

        const Wrapper = struct {
            fn render(ptr: *anyopaque, renderer: void) void {
                const self: Ptr = @ptrCast(@alignCast(ptr));
                self.render(renderer);
            }

            fn addChild(ptr: *anyopaque, child: Node) void {
                const self: Ptr = @ptrCast(@alignCast(ptr));
                self.addChild(child);
            }

            fn deinit(ptr: *anyopaque) void {
                const self: Ptr = @ptrCast(@alignCast(ptr));
                self.deinit();
            }
        };

        return .{
            .ptr = obj,
            .vtab = &.{
                .renderFn = Wrapper.render,
                .addChildFn = Wrapper.addChild,
                .deinitFn = Wrapper.deinit,
            },
            .type = obj.type,
        };
    }
};

// const Section = struct {
//     heading: Heading,
//     children: []Node(.Heading),
// };
// const Heading = struct {
//     level: u8,
//     title: []Node(.Text),
// };
// const Text = struct {
//     text: []u8,
//     children: ?[]Node(.Text | .Bold | .Italic),
// };
// const Bold = struct {
//     children: []Node(.Text | .Italic),
// };
// const Italic = struct {
//     children: []Node(.Text | .Bold),

//     pub fn addChild(self: *Italic, child: Node(.Text | .Bold)) void {
//         _ = self; // autofix
//         _ = child; // autofix
//         // TODO

//     }

//     pub fn toNode(self: *Italic) Node(.Italic) {
//         return Node.init(&self);
//     }
// };

pub const Document = struct {
    children: Array(Node),
    type: NodeType = .Document,
    pub fn toNode(self: *Document) Node {
        return Node.init(self);
    }

    var arena: std.heap.ArenaAllocator = undefined;
    var alloc: std.mem.Allocator = undefined;

    pub fn init(allocator: std.mem.Allocator) !*Document {
        arena = std.heap.ArenaAllocator.init(allocator);
        alloc = arena.allocator();
        const self = try alloc.create(Document);
        self.* = Document{
            .children = Array(Node).init(alloc),
        };
        return self;
    }
    pub fn deinit(self: *Document) void {
        arena.deinit();
        _ = self;
    }
    pub fn addChild(self: *Document, child: Node) void {
        _ = self; // autofix
        _ = child; // autofix
        // TODO
    }
    pub fn render(self: *Document, renderer: anytype) void {
        _ = self; // autofix
        _ = renderer; // autofix
        // TODO
    }
};

test "ASTNode init" {
    const alloc = std.testing.allocator;
    var arena = std.heap.ArenaAllocator.init(alloc);
    const a = arena.allocator();
    _ = a; // autofix
    var doc = try Document.init(alloc);
    var node = doc.toNode();
    defer node.deinit();

    std.debug.print("{any} \n", .{node});
}
