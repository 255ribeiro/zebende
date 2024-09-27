const std = @import("std");
const expect = std.testing.expect;
const pow = std.math.pow;
const Mutex = std.Thread.Mutex;
const RwLock = std.Thread.RwLock;
// c import

// const c = @cImport({
//     // @cInclude("stdio.h");
//     @cInclude("stdlib.h");
// });

// const testing = std.testing;

const Arr_writer = struct {
    value: [*c]f64,
    lock: RwLock = .{},
};

const print = std.debug.print;

pub const p_dcca_out = extern union {
    table: [*c][*c]f64,
    mat: [*c][*c][*c]f64,
};

const Output = extern struct {
    pdcca: p_dcca_out,
};

export fn p_dcca(
    data: [*c]f64,
    data_len: usize,
    data_count: usize,
    tws: [*c]usize,
    tws_len: usize,
    time_steps: [*c]f64,
    DCCA_of: [*c][*c]usize,
    DCCA_of_n_rows: usize,
    //outputs
    F_DFA_arr: [*c][*c]f64,
    DCCA_arr: [*c][*c]f64,
    P_DCCA_arr: [*c][*c]f64,
) void {

    // Allocation setup -- start
    print(" tipe of {}\n", .{@TypeOf(P_DCCA_arr)});

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Allocation setup -- end

    // Time scale loop -- start
    for (0..tws_len) |n_index| {
        const n = tws[n_index];
        const no_of_windows: usize = (data_len - n);

        // Allocation block

        var det_mat = allocator.alloc([][]f64, (n + 1)) catch |err| {
            print("Error allocating det_mat axis 0\n{}\n", .{err});
            break;
        };
        for (0..det_mat.len) |i| {
            det_mat[i] = allocator.alloc([]f64, no_of_windows) catch |err| {
                print("Error allocating det_mat axis 1 line{}\n{}\n", .{ i, err });
                break;
            };

            for (0..det_mat[0].len) |j| {
                det_mat[i][j] = allocator.alloc(f64, data_count) catch |err| {
                    print("Error allocating det_mat axis 2 line{} level {}\n{}\n", .{ i, j, err });
                    break;
                };
            }
        }

        defer {
            for (det_mat) |row| {
                for (row) |cel| {
                    allocator.free(cel);
                }
                allocator.free(row);
            }
            allocator.free(det_mat);
        }

        { // poll detrend scope -- start
            var pool: std.Thread.Pool = undefined;
            pool.init(.{ .allocator = allocator }) catch |err| {
                print("Error allocating thread pool {}\n", .{err});
                break;
            };
            defer pool.deinit();

            // windows loop -- start
            for (0..(no_of_windows * data_count)) |data_index| {
                pool.spawn(detend_prl, .{
                    &data,
                    &time_steps,
                    // usize vars
                    data_index,
                    no_of_windows,
                    n,
                    // mat output
                    &det_mat,
                }) catch |err| {
                    print("Error spawning thread pool no. {} error: {}\n", .{ data_index, err });
                    break;
                };
            } // windows loop -- end
        } // pool detrend scope -- end

        // dfa calc

        // { // pool F_DFA scope -- start
        //     var pool: std.Thread.Pool = undefined;
        //     pool.init(.{ .allocator = allocator }) catch |err| {
        //         print("Error allocating thread pool {}\n", .{err});
        //         break;
        //     };
        //     defer pool.deinit();

        for (0..(no_of_windows * data_count)) |data_index| {
            F_DFA_fill_prl(
                // input mat
                &det_mat,
                // usize vars
                data_index,
                n,
                n_index,
                data_count,
                // aoutput mat
                &F_DFA_arr,
            );
            // Parallel spawn -- start
            // pool.spawn(F_DFA_fill_prl, .{
            //     // input mat
            //     &det_mat,
            //     // usize vars
            //     data_index,
            //     n,
            //     n_index,
            //     data_count,
            //     // aoutput mat
            //     &F_DFA_arr,
            // }) catch |err| {
            //     print("Error spawning thread pool no. {} error: {}\n", .{ data_index, err });
            //     break;
            // };

            // Parallel spawn -- end
        }
        // } // pool F_DFA scope -- end

        for (0..data_count) |sr_ind| {
            F_DFA_arr[n_index][sr_ind] = @sqrt(F_DFA_arr[n_index][sr_ind] / @as(f64, @floatFromInt(no_of_windows)));
        }

        // dcca calc

        { // poll DCCA scope -- start
            var pool: std.Thread.Pool = undefined;
            pool.init(.{ .allocator = allocator }) catch |err| {
                print("Error allocating thread pool {}\n", .{err});
                break;
            };
            defer pool.deinit();

            for (0..(no_of_windows * DCCA_of_n_rows)) |data_index| {
                pool.spawn(DCCA_fill_prl, .{
                    // input matrices
                    &det_mat,   &DCCA_of,
                    // usize vars
                    data_index, n,
                    n_index,    DCCA_of_n_rows,
                    // aoutput mat
                    &DCCA_arr,
                }) catch |err| {
                    print("Error spawning thread pool no. {} error: {}\n", .{ data_index, err });
                    break;
                };
            }
        } // poll DCCA scope -- end

        for (0..DCCA_of_n_rows) |dcca_pair| {
            DCCA_arr[n_index][dcca_pair] = DCCA_arr[n_index][dcca_pair] / @as(f64, @floatFromInt(no_of_windows));
        }
        // p_dcca calc

        for (0..DCCA_of_n_rows) |dcca_pair| {
            DCCA_arr[n_index][dcca_pair] = DCCA_arr[n_index][dcca_pair] / @as(f64, @floatFromInt(no_of_windows));

            p_dcca_simple_output(
            // mat imputs
            &F_DFA_arr, &DCCA_arr, &DCCA_of,
            // usize paris
            dcca_pair, n_index,
            // output
            &P_DCCA_arr);
        }

        // Allocation free
    } // Time scale loop -- end
} // function --end

pub fn lin_ls_fit(win: []f64, time: []f64) [2]f64 {
    var x_sum: f64 = 0;
    var y_sum: f64 = 0;
    var xy_sum: f64 = 0;
    var x2_sum: f64 = 0;
    for (win, time) |w, t| {
        x_sum += t;
        y_sum += w;
        xy_sum += t * w;
        x2_sum += pow(f64, t, 2);
    }
    const n: f64 = @as(f64, @floatFromInt(time.len));
    //slope
    const slope: f64 = (((n * xy_sum) - (x_sum * y_sum)) /
        ((n * x2_sum) - (pow(f64, x_sum, 2))));
    //inter
    const inter: f64 = ((y_sum - (slope * x_sum)) /
        (n));
    //result
    return [_]f64{ slope, inter };
}

pub fn detend_prl(
    // data input
    data: *const [*c]f64,
    time_steps: *const [*c]f64,
    // usize vars
    data_index: usize,
    no_of_windows: usize,
    n: usize,
    // mat output
    det_mat: *const [][][]f64,
) void {
    const tw_start: usize = data_index % no_of_windows;
    const sr_ind = @divFloor(data_index, no_of_windows);
    const sw_start: usize = data_index + (sr_ind * n);

    const time = time_steps.*[tw_start .. tw_start + (n + 1)];
    const win = data.*[sw_start .. sw_start + (n + 1)];
    const slp_int = lin_ls_fit(win, time);
    for (win, time, 0..) |w, t, i| {
        det_mat.*[i][tw_start][sr_ind] = w - (slp_int[0] * t + slp_int[1]);
    }
}

pub fn F_DFA_fill_prl(
    // input mat
    det_mat: *const [][][]f64,
    // usize vars
    data_index: usize,
    n: usize,
    n_index: usize,
    data_count: usize,
    // aoutput mat
    F_DFA_arr: *const [*c][*c]f64,
) void {
    const sr_index = data_index % data_count;
    const win_index = @divFloor(data_index, data_count);

    var arr_writer = Arr_writer{ .value = &F_DFA_arr.*[n_index][sr_index], .lock = .{} };

    var dfa_temp: f64 = 0;
    for (0..(n + 1)) |i| {
        dfa_temp += pow(f64, det_mat.*[i][win_index][sr_index], 2);
    }

    arr_writer.lock.lock();
    defer arr_writer.lock.unlock();
    arr_writer.value.* += dfa_temp / @as(f64, @floatFromInt(n + 1));
}

pub fn DCCA_fill_prl(
    // input matrices
    det_mat: *const [][][]f64,
    DCCA_of: *const [*c][*c]usize,
    // usize vars
    data_index: usize,
    n: usize,
    n_index: usize,
    DCCA_of_n_rows: usize,
    // aoutput mat
    DCCA_arr: *const [*c][*c]f64,
) void {
    const dcca_pair = data_index % DCCA_of_n_rows;
    const win_index = @divFloor(data_index, DCCA_of_n_rows);
    var dcca_temp: f64 = 0;
    for (0..(n + 1)) |i| {
        dcca_temp += det_mat.*[i][win_index][DCCA_of.*[dcca_pair][0]] * det_mat.*[i][win_index][DCCA_of.*[dcca_pair][1]];
    }
    DCCA_arr.*[n_index][dcca_pair] = dcca_temp / @as(f64, @floatFromInt(n + 1));
}

pub fn dcca_from_columns(
    det_mat: [][][]f64,
    DCCA_of: *const [*c][*c]usize,
    dcca_pair: usize,
    n_index: usize,

    // output
    DCCA_arr: *const [*c][*c]f64,
) void {
    const place_holder: usize = 0;
    var dcca_temp: f64 = 0;

    for (0..det_mat.len) |j| {
        dcca_temp += det_mat[j][DCCA_of.*[dcca_pair][0]][place_holder] * det_mat[j][DCCA_of.*[dcca_pair][1]][place_holder];
    }
    DCCA_arr.*[n_index][dcca_pair] += dcca_temp / @as(f64, @floatFromInt(det_mat.len));
}

pub fn p_dcca_simple_output(
    // mat imputs
    F_DFA_arr: *const [*c][*c]f64,
    DCCA_arr: *const [*c][*c]f64,
    // dcca of
    DCCA_of: *const [*c][*c]usize,
    // usize vars
    dcca_pair: usize,
    n_index: usize,
    // output
    P_DCCA_arr: *const [*c][*c]f64,
) void {
    P_DCCA_arr.*[n_index][dcca_pair] = DCCA_arr.*[n_index][dcca_pair] / (F_DFA_arr.*[n_index][DCCA_of.*[dcca_pair][0]] * F_DFA_arr.*[n_index][DCCA_of.*[dcca_pair][1]]);
}

pub fn p_dcca_matrix_output(
    // mat imputs
    F_DFA_arr: *const [*c][*c]f64,
    DCCA_arr: *const [*c][*c]f64,
    // dcca of
    DCCA_of: *const [*c][*c]usize,
    // usize vars
    dcca_pair: usize,
    n_index: usize,
    // output
    P_DCCA_arr: *const [*c][*c][*c]f64,
) void {
    print("{}, {}", .{ n_index, F_DFA_arr, DCCA_arr, DCCA_of, dcca_pair, P_DCCA_arr });
}
