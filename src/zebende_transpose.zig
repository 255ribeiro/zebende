const std = @import("std");
const expect = std.testing.expect;
const pow = std.math.pow;
const Now_ms = std.time.milliTimestamp;
// c import

// const c = @cImport({
//     // @cInclude("stdio.h");
//     @cInclude("stdlib.h");
// });

// const testing = std.testing;

const print = std.debug.print;

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
) void {
    print(" tipe of {}\n", .{@TypeOf(P_DCCA_arr)});

    // Allocation setup -- start
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    // Allocation setup -- end

    // Allocation block -- start
    var arena = std.heap.ArenaAllocator.init(gpa.allocator());
    const allocator = arena.allocator();
    defer arena.deinit();

    var dcca_n = allocator.alloc([]f64, DCCA_of_n_rows) catch unreachable;

    var f2dfa_n = allocator.alloc([]f64, data_count) catch unreachable;

    var det_mat = allocator.alloc([]f64, data_count) catch unreachable;

    for (0..dcca_n.len) |i| {
        dcca_n[i] = allocator.alloc(f64, (data_len - tws[0])) catch |err| {
            print("Error allocating dcca_n axis 0\n{}\n", .{err});
            break;
        };
    }

    for (0..f2dfa_n.len) |i| {
        f2dfa_n[i] = allocator.alloc(f64, (data_len - tws[0])) catch |err| {
            print("Error allocating f2dfa_n axis 1 line{}\n{}\n", .{ i, err });
            break;
        };
    }

    for (0..det_mat.len) |i| {
        det_mat[i] = allocator.alloc(f64, (tws[tws_len - 1] + 1)) catch |err| {
            print("Error allocating det_mat axis 1 line{}\n{}\n", .{ i, err });
            break;
        };
    }

    // Allocation block -- end

    // Time scale loop -- sart
    for (0..tws_len, tws[0..]) |n_index, n| {
        const no_of_windows = (data_len - n);
        const window_size = n + 1;
        // for each Time window
        // windows loop -- start
        for (0..data_len - n) |win_start| {
            // print("wind_size: {} - window start: {}\n", .{ n, win_start });

            // time windows slice
            const time_window = time_steps[win_start .. win_start + (n + 1)];

            // for each data series
            // series loop -- start
            for (0..data_count) |sr_ind| {
                const serie = data[(sr_ind * data_len)..((sr_ind * data_len) + data_len)];
                // print("---\nrow num: {d} row size: {d} row: {d}\n---\n", .{ sr_ind, serie.len, serie });

                // data slice
                const window = serie[win_start .. win_start + (n + 1)];

                // print("window: {d}\nTime window: {d}\n", .{ window, time_window });

                detrend_f2dfa_n(window, time_window, sr_ind, win_start, det_mat, f2dfa_n);
            } // series loop -- end

            for (0..DCCA_of_n_rows) |dcca_pair| {
                dcca_from_columns(det_mat, DCCA_of, dcca_pair, window_size, win_start, dcca_n);
            }
        } // windows loop -- end

        for (0..data_count) |sr_ind| {
            // Fill dfa table
            fill_F_DFA_arr(
                n_index,
                f2dfa_n,
                sr_ind,
                no_of_windows,
                // output
                F_DFA_arr,
            );
        }

        for (0..DCCA_of_n_rows) |dcca_pair| {
            // fill DCCA  table
            fill_DCCA_arr(
                n_index,
                dcca_n,
                no_of_windows,
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

        // Allocation free

    } // time window scale -- loop end

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
    sr_ind: usize,
    win_start: usize,
    det_mat: [][]f64,
    f2dfa_n: [][]f64,
) void {
    const slp_int = linLestSquaresFit(win, time);
    var dfa_temp: f64 = 0;
    for (win, time, 0..) |w, t, i| {
        det_mat[sr_ind][i] = w - (slp_int[0] * t + slp_int[1]);
        dfa_temp += pow(f64, det_mat[sr_ind][i], 2);
    }
    f2dfa_n[sr_ind][win_start] = dfa_temp / @as(f64, @floatFromInt(time.len));
}

pub fn dcca_from_columns(
    det_mat: [][]f64,
    DCCA_of: [*c][*c]usize,
    dcca_pair: usize,
    window_size: usize,
    win_start: usize,
    dcca_n: [][]f64,
) void {
    var dcca_temp: f64 = 0;
    for (0..window_size) |j| {
        dcca_temp += det_mat[DCCA_of[dcca_pair][0]][j] * det_mat[DCCA_of[dcca_pair][1]][j];
    }
    dcca_n[dcca_pair][win_start] = dcca_temp / @as(f64, @floatFromInt(window_size));
}

pub fn fill_F_DFA_arr(n_index: usize, f2dfa_n: [][]f64, sr_ind: usize, no_of_windows: usize, F_DFA_arr: [*c][*c]f64) void {
    var dfa_temp: f64 = 0;

    for (0..no_of_windows) |win| {
        dfa_temp += f2dfa_n[sr_ind][win];
    }

    F_DFA_arr[n_index][sr_ind] = @sqrt(dfa_temp / @as(f64, @floatFromInt(no_of_windows)));
}

pub fn fill_DCCA_arr(
    n_index: usize,
    dcca_n: [][]f64,
    no_of_windows: usize,
    dcca_pair: usize,
    // output
    DCCA_arr: [*c][*c]f64,
) void {
    var dcca_temp: f64 = 0;

    for (0..no_of_windows) |win| {
        dcca_temp += dcca_n[dcca_pair][win];
    }

    DCCA_arr[n_index][dcca_pair] = dcca_temp / @as(f64, @floatFromInt(no_of_windows));
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
