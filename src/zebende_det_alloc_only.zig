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

    // Time scale loop -- start
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

        const start_det_time = Now_ms();
        const no_of_windows: usize = (data_len - n);
        // windows loop -- start
        for (0..no_of_windows) |window_start| {
            // time windows slice
            const time_window = time_steps[window_start .. window_start + (n + 1)];

            // for each data series
            // series loop -- start
            for (0..data_count) |sr_ind| {
                const serie = data[(sr_ind * data_len)..((sr_ind * data_len) + data_len)];

                // data slice
                const window = serie[window_start .. window_start + (n + 1)];

                detrend_dfa(
                // slice inputs
                window, time_window,
                // usize inputs
                sr_ind, n_index,

                // outputs
                det_mat, &F_DFA_arr);
            } // series loop -- end

            for (0..DCCA_of_n_rows) |dcca_pair| {
                dcca_from_columns(
                // mat inputs
                det_mat, &DCCA_of,
                // usize var inputs
                dcca_pair, n_index,
                // output
                &DCCA_arr);
            }
        } // windows loop -- end

        const elapsed_det_time = Now_ms() - start_det_time;
        print("elapsed det time {}\n", .{elapsed_det_time});
        //
        const start_filldfa_time = Now_ms();

        for (0..data_count) |sr_ind| {
            F_DFA_arr[n_index][sr_ind] = @sqrt(F_DFA_arr[n_index][sr_ind] / @as(f64, @floatFromInt(no_of_windows)));
        }

        const elapsed_filldfa_time = Now_ms() - start_filldfa_time;
        print("elapsed fill dfa time {}\n", .{elapsed_filldfa_time});

        const start_dcca_pdcca_time = Now_ms();

        for (0..DCCA_of_n_rows) |dcca_pair| {
            DCCA_arr[n_index][dcca_pair] = DCCA_arr[n_index][dcca_pair] / @as(f64, @floatFromInt(no_of_windows));

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

        const tws_elapsed = (Now_ms() - start_t_loop);
        print("tws elapsed: {} milis for time scale no. {} of size {}\n", .{ tws_elapsed, n_index, n });
        // Allocation free
    } // Time scale loop -- end
    print("elapsed time for all {}\n", .{Now_ms() - start_time});
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

pub fn detrend_dfa(
    win: []f64,
    time: []f64,
    sr_ind: usize,
    n_index: usize,
    // output
    det_mat: [][]f64,
    F_DFA_arr: *const [*c][*c]f64,
) void {
    const slp_int = lin_ls_fit(win, time);
    var dfa_temp: f64 = 0;
    for (win, time, 0..) |w, t, i| {
        det_mat[i][sr_ind] = w - (slp_int[0] * t + slp_int[1]);
        dfa_temp += pow(f64, det_mat[i][sr_ind], 2);
    }
    F_DFA_arr.*[n_index][sr_ind] += dfa_temp / @as(f64, @floatFromInt(time.len));
}

pub fn dcca_from_columns(
    det_mat: [][]f64,
    DCCA_of: *const [*c][*c]usize,
    dcca_pair: usize,
    n_index: usize,

    // output
    DCCA_arr: *const [*c][*c]f64,
) void {
    var dcca_temp: f64 = 0;
    for (0..det_mat.len) |j| {
        dcca_temp += det_mat[j][DCCA_of.*[dcca_pair][0]] * det_mat[j][DCCA_of.*[dcca_pair][1]];
    }
    DCCA_arr.*[n_index][dcca_pair] += dcca_temp / @as(f64, @floatFromInt(det_mat.len));
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
