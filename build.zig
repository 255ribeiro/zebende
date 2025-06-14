const std = @import("std");

const lib_name = "zebendezig";

const targets: []const std.Target.Query = &.{
    .{ .cpu_arch = .aarch64, .os_tag = .macos },
    .{ .cpu_arch = .aarch64, .os_tag = .linux },
    .{ .cpu_arch = .aarch64, .os_tag = .windows },
    .{ .cpu_arch = .x86_64, .os_tag = .windows },
    .{ .cpu_arch = .x86_64, .os_tag = .macos },
    .{ .cpu_arch = .x86_64, .os_tag = .linux },
};

pub fn build(b: *std.Build) !void {
    var src_file_name: []const u8 = "zebende_opt";

    const ball = b.option(bool, "ball", "build for all platforms") orelse false;
    const basic = b.option(bool, "basic", "build basic algorithm") orelse false;

    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});
    if (basic) {
        src_file_name = "zebende_basic";
    }

    // creating path string
    const allocator = std.heap.page_allocator;
    const src_path = try std.mem.concat(allocator, u8, &.{ "src/", src_file_name, ".zig" });
    defer allocator.free(src_path);

    if (ball) {
        for (targets) |t| {
            const lib = b.addSharedLibrary(.{
                .name = lib_name,
                .root_source_file = b.path(src_path),
                .target = b.resolveTargetQuery(t),
                .optimize = optimize,
            });

            const target_output = b.addInstallArtifact(lib, .{
                .dest_dir = .{
                    .override = .{
                        .custom = try t.zigTriple(b.allocator),
                    },
                },
            });

            b.getInstallStep().dependOn(&target_output.step);
        }
    } else {
        const lib = b.addSharedLibrary(.{
            .target = target,
            .optimize = optimize,
            .name = lib_name,
            .root_source_file = b.path(src_path),
        });

        const target_output = b.addInstallArtifact(lib, .{
            .dest_dir = .{
                .override = .{
                    .custom = "./x86_64-windows/",
                },
            },
        });

        b.getInstallStep().dependOn(&target_output.step);
    }
    std.debug.print("compiling: {s}\n", .{src_path});
}
