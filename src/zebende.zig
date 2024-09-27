const std = @import("std");
const print = std.debug.print;
const expect = std.testing.expect;
const pow = std.math.pow;
const Now_ms = std.time.milliTimestamp;
// const testing = std.testing;

// c import
// const c = @cImport({
//     // @cInclude("stdio.h");
//     @cInclude("stdlib.h");
// });

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
    P_DCCA_arr: [*c]f64,
    P_dcca_matrix_output_flag: bool,
    // P_dcca_matrix_flag: bool,
) void {
    const start_time = Now_ms();
    print("Starting calculations:\n\nNumber of time steps:\t\t{}\nNumber of series:\t\t{}\nSize of series:\t\t\t{}\nNumber of DCCA combinations:\t{}\n\n", .{
        // initial values
        tws_len, data_count, data_len, DCCA_of_n_rows,
    });

    // Allocation setup -- start
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    // Allocation setup -- end

    // Time scale loop -- sart
    for (0..tws_len, tws[0..]) |n_index, n| {
        const start_t_loop = Now_ms();
        print("\n---> starting calculations for time scale no. {} of size {}\n", .{
            n_index,
            n,
        });
        const start_alloc_time = Now_ms();
        // Allocation block -- start
        var arena = std.heap.ArenaAllocator.init(gpa.allocator());
        const allocator = arena.allocator();
        defer arena.deinit();

        var dcca_n = allocator.alloc([]f64, (data_len - n)) catch |err| {
            print("Error allocating dcca_n axis 0\n{}\n", .{err});
            break;
        };
        for (0..dcca_n.len) |i| {
            dcca_n[i] = allocator.alloc(f64, DCCA_of_n_rows) catch |err| {
                print("Error allocating dcca_n axis 0\n{}\n", .{err});
                break;
            };
        }

        var f2dfa_n = allocator.alloc([]f64, (data_len - n)) catch |err| {
            print("Error allocating f2dfa_n axis 0\n{}\n", .{err});
            break;
        };
        for (0..f2dfa_n.len) |i| {
            f2dfa_n[i] = allocator.alloc(f64, data_count) catch |err| {
                print("Error allocating f2dfa_n axis 1 line{}\n{}\n", .{ i, err });
                break;
            };
        }

        var det_mat = allocator.alloc([]f64, (n + 1)) catch |err| {
            print("Error allocating det_mat axis 0\n{}\n", .{err});
            break;
        };
        for (0..det_mat.len) |i| {
            det_mat[i] = allocator.alloc(f64, data_count) catch |err| {
                print("Error allocating det_mat axis 1 line{}\n{}\n", .{ i, err });
                break;
            };
        }
        // Allocation block -- end
        const elapsed_alloc_time = Now_ms() - start_alloc_time;
        print("elapsed allocation time {}\n", .{elapsed_alloc_time});

        // for each Time window
        // time windows loop -- start

        const start_det_time = Now_ms();
        for (0..data_len - n) |win_start| {
            // print("wind_size: {} - window start: {}\n", .{ n, win_start });

            // time windows slice
            const time_window = time_steps[win_start..][0..(n + 1)];

            // for each data series
            // series loop -- start
            for (0..data_count) |sr_index| {
                const serie = data[(sr_index * data_len)..][0..data_len];
                // print("---\nrow num: {d} row size: {d} row: {d}\n---\n", .{ sr_ind, serie.len, serie });

                // data slice
                const window = serie[win_start..][0..(n + 1)];

                // print("window: {d}\nTime window: {d}\n", .{ window, time_window });

                detrend_f2dfa_n(window, time_window, sr_index, win_start, det_mat, f2dfa_n);
            } // series loop -- end

            for (0..DCCA_of_n_rows) |dcca_pair| {
                dcca(det_mat, DCCA_of, dcca_pair, win_start, dcca_n);
            }
        } // time windows loop -- end
        const elapsed_det_time = Now_ms() - start_det_time;
        print("elapsed det time {}\n", .{elapsed_det_time});
        //
        const start_filldfa_time = Now_ms();
        for (0..data_count) |sr_index| {
            // Fill dfa table
            fill_F_DFA_arr(
                n_index,
                f2dfa_n,
                sr_index,
                // output
                F_DFA_arr,
            );
        }

        const elapsed_filldfa_time = Now_ms() - start_filldfa_time;
        print("elapsed fill dfa time {}\n", .{elapsed_filldfa_time});

        const start_dcca_pdcca_time = Now_ms();
        for (0..DCCA_of_n_rows) |dcca_pair| {
            // fill DCCA  table
            fill_DCCA_arr(
                n_index,
                dcca_n,
                dcca_pair,
                // output
                DCCA_arr,
            );
            if (P_dcca_matrix_output_flag == false) {
                pdccaSimpleOutput(n_index, F_DFA_arr, DCCA_arr, DCCA_of, DCCA_of_n_rows, dcca_pair,
                // output
                P_DCCA_arr);
            } else {
                pdccaMatrixOutput(n_index, tws_len, data_count, F_DFA_arr, DCCA_arr, DCCA_of, dcca_pair,
                // output
                P_DCCA_arr);
            }
        }
        const elapsed_pdcca_time = Now_ms() - start_dcca_pdcca_time;
        print("elapsed DCCA P_DCCA time {}\n", .{elapsed_pdcca_time});
        // Allocation free
        const tws_elapsed = (Now_ms() - start_t_loop);
        print("tws elapsed: {} milis for time scale no. {} of size {}\n", .{ tws_elapsed, n_index, n });
    } // time window scale -- loop end

    print("elapsed time for all {}\n", .{Now_ms() - start_time});
} // function -- end

pub fn linLestSquaresFit(win: []f64, time: []f64) [2]f64 {
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

pub fn detrend_f2dfa_n(
    win: []f64,
    time: []f64,
    sr_index: usize,
    win_start: usize,
    det_mat: [][]f64,
    f2dfa_n: [][]f64,
) void {
    const slp_int = linLestSquaresFit(win, time);
    var dfa_temp: f64 = 0;
    for (win, time, 0..) |w, t, i| {
        det_mat[i][sr_index] = w - (slp_int[0] * t + slp_int[1]);
        dfa_temp += pow(f64, det_mat[i][sr_index], 2);
    }
    f2dfa_n[win_start][sr_index] = dfa_temp / @as(f64, @floatFromInt(time.len));
}

pub fn dcca(
    det_mat: [][]f64,
    DCCA_of: [*c][*c]usize,
    dcca_pair: usize,
    win_start: usize,
    dcca_n: [][]f64,
) void {
    var dcca_temp: f64 = 0;
    for (0..det_mat.len) |j| {
        dcca_temp += det_mat[j][DCCA_of[dcca_pair][0]] * det_mat[j][DCCA_of[dcca_pair][1]];
    }
    dcca_n[win_start][dcca_pair] = dcca_temp / @as(f64, @floatFromInt(det_mat.len));
}

pub fn fill_F_DFA_arr(n_index: usize, f2dfa_n: [][]f64, sr_index: usize, F_DFA_arr: [*c][*c]f64) void {
    var dfa_temp: f64 = 0;

    for (0..f2dfa_n.len) |win| {
        dfa_temp += f2dfa_n[win][sr_index];
    }

    F_DFA_arr[n_index][sr_index] = @sqrt(dfa_temp / @as(f64, @floatFromInt(f2dfa_n.len)));
}

pub fn fill_DCCA_arr(
    n_index: usize,
    dcca_n: [][]f64,
    dcca_pair: usize,
    // output
    DCCA_arr: [*c][*c]f64,
) void {
    var dcca_temp: f64 = 0;

    for (0..dcca_n.len) |win| {
        dcca_temp += dcca_n[win][dcca_pair];
    }

    DCCA_arr[n_index][dcca_pair] = dcca_temp / @as(f64, @floatFromInt(dcca_n.len));
}

pub fn pdccaSimpleOutput(
    n_index: usize,
    F_DFA_arr: [*c][*c]f64,
    DCCA_arr: [*c][*c]f64,
    // dcca of
    DCCA_of: [*c][*c]usize,
    DCCA_of_n_rows: usize,
    dcca_pair: usize,
    // output
    P_DCCA_arr: [*c]f64,
) void {
    P_DCCA_arr[n_index * DCCA_of_n_rows + dcca_pair] = DCCA_arr[n_index][dcca_pair] / (F_DFA_arr[n_index][DCCA_of[dcca_pair][0]] * F_DFA_arr[n_index][DCCA_of[dcca_pair][1]]);
}

pub fn pdccaMatrixOutput(
    n_index: usize,
    tws_len: usize,
    data_count: usize,
    F_DFA_arr: [*c][*c]f64,
    DCCA_arr: [*c][*c]f64,

    // dcca of
    DCCA_of: [*c][*c]usize,
    dcca_pair: usize,
    // output
    P_DCCA_arr: [*c]f64,
) void {

    // index = i * (dim2 * dim3) + j * dim3 + k;
    const pdcca = DCCA_arr[n_index][dcca_pair] / (F_DFA_arr[n_index][DCCA_of[dcca_pair][0]] * F_DFA_arr[n_index][DCCA_of[dcca_pair][1]]);

    P_DCCA_arr[DCCA_of[dcca_pair][0] * (data_count * tws_len) + DCCA_of[dcca_pair][1] * tws_len + n_index] = pdcca;
    P_DCCA_arr[DCCA_of[dcca_pair][1] * (data_count * tws_len) + DCCA_of[dcca_pair][0] * tws_len + n_index] = pdcca;
}
