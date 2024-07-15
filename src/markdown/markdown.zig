const std = @import("std");
const AST = @import("ast.zig");

const TokenType = enum {
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

const Token = struct {
    type: TokenType,
    loc: usize,
    ptr: [*]const u8,
    len: usize = 0,
};

const Needles = std.ComptimeStringMap(TokenType, .{
    .{ "#", .Heading },
    .{ "**", .Bold },
    .{ "_", .Italic },
    .{ "`", .Code },
    .{ "- ", .ListItem },
    .{ "> ", .BlockQuote },
    .{ "---", .HorizontalRule },
});

pub fn tokenize(alloc: std.mem.Allocator, input: []const u8) ![]Token {
    var tokens = std.ArrayList(Token).init(alloc);
    defer tokens.deinit();
    const startsWith = std.mem.startsWith;
    var idx: usize = 0;
    while (idx < input.len) : (idx += 1) {
        inline for (Needles.kvs) |kv| {
            const needle = kv.key;
            const tokenType = kv.value;
            if (startsWith(u8, input[idx..], needle)) {
                const token = Token{ .type = tokenType, .loc = idx, .ptr = input[idx + 1 ..].ptr };
                try tokens.append(token);
                idx += needle.len;
                break;
            }
        }
    }

    return tokens.toOwnedSlice();
}

pub fn parseToAST(alloc: std.mem.Allocator, tokens: []Token) !AST.Node {
    _ = tokens; // autofix
    var doc = try AST.Document.init(alloc);
    return doc.toNode();
}

fn parseToANSI(alloc: std.mem.Allocator, tokens: []Token) ![]u8 {
    var output = std.ArrayList(u8).init(alloc);
    defer output.deinit();

    for (tokens) |tok| {
        switch (tok.type) {
            TokenType.Heading => {
                try output.appendSlice("\x1b[1m");
            },
            TokenType.Text => {
                try output.appendSlice(tok.ptr[0..tok.len]);
            },
            TokenType.Bold => {
                try output.appendSlice("\x1b[1m");
            },
            TokenType.Italic => {
                try output.appendSlice("\x1b[3m");
            },
            TokenType.Code => {
                try output.appendSlice("\x1b[2m");
            },
            TokenType.ListItem => {
                try output.appendSlice("  • ");
            },
            TokenType.Paragraph => {
                try output.appendSlice("\n\n");
            },
            TokenType.Link => {
                try output.appendSlice("\x1b[4m");
            },
            TokenType.Image => {
                try output.appendSlice("\x1b[5m");
            },
            TokenType.BlockQuote => {
                try output.appendSlice("  > ");
            },
            TokenType.HorizontalRule => {
                try output.appendSlice("\n\n");
            },
        }
    }

    return output.toOwnedSlice();
}

const testText =
    \\# Title
    \\
    \\## Subtitle
    \\
    \\This is a paragraph with some text. It has a [link](https://example.com) and an image: ![alt text](https://example.com/image.png)
    \\
    \\- This is a list item
    \\- This is another list item
    \\
    \\---
    \\> This is a block quote
    \\
    \\**bold text**
    \\_italic text_
    \\`code text`
    \\```code block
    \\some code
    \\```
    \\
;

test "tokenize" {
    const allocator = std.testing.allocator;
    const tokens = try tokenize(allocator, testText);
    defer allocator.free(tokens);
    // std.debug.print("Tokens: {any}\n", .{tokens});

    var doc = try parseToAST(allocator, tokens);
    defer doc.deinit();
    std.debug.print("Doc: {any}\n", .{doc});
}
