const std = @import("std");

// const src_file_name = "zebende_basic";
const src_file_name = "zebende_opt";

const lib_name = "zebendezig";

const targets: []const std.Target.Query = &.{
    .{ .cpu_arch = .aarch64, .os_tag = .macos },
    .{ .cpu_arch = .aarch64, .os_tag = .linux },
    .{ .cpu_arch = .aarch64, .os_tag = .windows },
    .{ .cpu_arch = .x86_64, .os_tag = .windows },
    .{ .cpu_arch = .x86_64, .os_tag = .macos },
    .{ .cpu_arch = .x86_64, .os_tag = .linux },
    // .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .gnu },
    // .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .musl },
};

pub fn build(b: *std.Build) !void {
    const ball = b.option(bool, "ball", "build for all platforms") orelse false;
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    if (ball) {
        for (targets) |t| {
            const lib = b.addSharedLibrary(.{
                .name = lib_name,
                .root_source_file = b.path("src/" ++ src_file_name ++ ".zig"),
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
            .root_source_file = b.path("src/" ++ src_file_name ++ ".zig"),
        });

        b.installArtifact(lib);
    }
}
