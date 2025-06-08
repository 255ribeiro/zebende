const std = @import("std");

// const src_file_name = "zebende_basic";
const src_file_name = "zebende_opt";

const targets: []const std.Target.Query = &.{
    .{ .cpu_arch = .aarch64, .os_tag = .macos },
    .{ .cpu_arch = .aarch64, .os_tag = .linux },
    .{ .cpu_arch = .aarch64, .os_tag = .windows },
    .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .gnu },
    .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .musl },
    .{ .cpu_arch = .x86_64, .os_tag = .windows },
    .{ .cpu_arch = .x86_64, .os_tag = .macos },
};

pub fn build(b: *std.Build) !void {
    for (targets) |t| {
        const lib = b.addSharedLibrary(.{
            .name = src_file_name,
            .root_source_file = b.path("src/" ++ src_file_name ++ ".zig"),
            .target = b.resolveTargetQuery(t),
            .optimize = .ReleaseSafe,
        });
        // lib.setOutputDir("./zebende-out");
        const target_output = b.addInstallArtifact(lib, .{
            .dest_dir = .{
                .override = .{
                    .custom = try t.zigTriple(b.allocator),
                },
            },
        });

        b.getInstallStep().dependOn(&target_output.step);
    }
}
