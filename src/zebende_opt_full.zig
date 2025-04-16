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

const TimeData = struct {
    sum_x: f64 = 0,
    sum_x2: f64 = 0,
    window_len: usize = 0,

    fn init(allocator: std.mem.Allocator) *TimeData {
        const ptr = allocator.create(TimeData) catch unreachable;
        ptr.sum_x = 0;
        ptr.sum_x2 = 0;
        ptr.window_len = 0;
        return ptr;
    }
};

const DCCAOperator = struct {
    previous: *TimeData,
    current: *TimeData,
    time: []f64,
    left_x: f64,
    time_window: []f64,
    series: []*SerieContainer,
    allocator: std.mem.Allocator,

    fn init(allocator: std.mem.Allocator, time: []f64, series: *[]*SerieContainer) *DCCAOperator {
        const ptr = allocator.create(DCCAOperator) catch unreachable;
        ptr.previous = TimeData.init(allocator);
        ptr.current = TimeData.init(allocator);
        ptr.time = time;
        allocator = allocator;

        ptr.series = series;

        return ptr;
    }
    // shift window
    fn shiftWindow(self: *DCCAOperator, n: usize, win_start: usize, F_DFA_ptr: *allowzero [*c]f64) void {
        self.time_window = self.time[win_start .. win_start + (n + 1)];
        // print("win_start {}\n", .{win_start});
        if (win_start != 0) {
            self.current.sum_x = self.current.sum_x - self.left_x + self.time_window[n];

            self.current.sum_x2 = self.current.sum_x2 - pow(f64, self.left_x, 2) + pow(f64, self.time_window[n], 2);

            for (self.series, 0..) |serie, sr_index| {
                serie.current.sum_y = serie.current.sum_y - serie.left_y + serie.serie[win_start + n];
                serie.current.sum_xy = serie.current.sum_xy - (self.left_x * serie.left_y) + (self.time_window[n] * serie.serie[win_start + n]);

                serie.left_y = serie.serie[win_start];

                self.detrended(serie, win_start, &F_DFA_ptr.*[sr_index]);
            }
        } else { // win_start == 0
            self.current.window_len = n + 1;

            for (self.previous.window_len..self.current.window_len) |i| {
                self.previous.sum_x += self.time_window[i];
                self.previous.sum_x2 += pow(f64, self.time_window[i], 2);

                for (self.series) |serie| {
                    serie.previous.sum_y += serie.serie[i];
                    serie.previous.sum_xy += self.time_window[i] * serie.serie[i];
                }
            }

            // updating current sum values
            self.current.sum_x = self.previous.sum_x;
            self.current.sum_x2 = self.previous.sum_x2;

            for (self.series, 0..) |serie, sr_index| {
                serie.current.sum_xy = serie.previous.sum_xy;
                serie.current.sum_y = serie.previous.sum_y;

                serie.left_y = serie.serie[win_start];
                self.detrended(serie, win_start, &F_DFA_ptr.*[sr_index]);
            }
        }
        self.left_x = self.time_window[0];
    }

    //deternded function

    fn detrended(self: *DCCAOperator, serie: *SerieContainer, win_start: usize, F_DFA_ptr: *allowzero f64) void {
        const len: f64 = @as(f64, @floatFromInt(self.current.window_len));
        //slope
        serie.slope = (((len * serie.current.sum_xy) - (self.current.sum_x * serie.current.sum_y)) /
            ((len * self.current.sum_x2) - (pow(f64, self.current.sum_x, 2))));
        //inter
        serie.inter = ((serie.current.sum_y - (serie.slope * self.current.sum_x)) /
            (len));
        //result
        serie.dfa_temp = 0;
        for (self.time_window, 0..) |t, i| {
            serie.det_row.*[i] = serie.serie[win_start + i] - (serie.slope * t + serie.inter);
            serie.dfa_temp += pow(f64, serie.det_row.*[i], 2);
        }
        F_DFA_ptr.* += serie.dfa_temp / len;
    }

    // reset function
    fn reset(
        self: DCCAOperator,
    ) void {
        self.previous.window_len = self.current.window_len;
    }
};

const SerieData = struct {
    sum_y: f64 = 0,
    sum_xy: f64 = 0,

    fn init(allocator: std.mem.Allocator) *SerieData {
        const ptr = allocator.create(SerieData) catch unreachable;
        ptr.sum_xy = 0;
        ptr.sum_y = 0;

        return ptr;
    }
};

const SerieContainer = struct {
    previous: *SerieData,
    current: *SerieData,
    serie: []f64,
    det_row: *const []f64,
    left_y: f64,
    slope: f64,
    inter: f64,
    dfa_temp: f64,

    fn init(allocator: std.mem.Allocator) *SerieContainer {
        const ptr = allocator.create(SerieContainer) catch unreachable;
        ptr.previous = SerieData.init(allocator);
        ptr.current = SerieData.init(allocator);
        return ptr;
    }
};

export fn p_dcca(
    data: [*c]f64,
    data_len: usize,
    data_count: usize,
    tws: [*c]usize,
    tws_len: usize,
    time_steps: [*c]f64,
    //outputs
    F_DFA_arr: [*c][*c]f64,
    DCCA_arr: [*c][*c]f64,
    P_DCCA_arr: [*c]f64,
    // p_dcca matrix output flag
    P_dcca_matrix_output_flag: bool,
) void {
    const start_time = Now_ms();
    print("Starting calculations:\n\nNumber of time steps:\t\t{}\nNumber of series:\t\t{}\nSize of series:\t\t\t{}\nNumber of DCCA combinations:\t{}\n\n", .{
        // initial values
        tws_len, data_count, data_len,
    });

    // Allocation setup -- start
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var arena = std.heap.ArenaAllocator.init(gpa.allocator());
    const allocator = arena.allocator();
    defer arena.deinit();
    // Allocation setup -- end

    var det_mat = allocator.alloc([]f64, data_count) catch unreachable;

    var series_array = allocator.alloc(*SerieContainer, data_count) catch unreachable;

    for (0..data_count) |sr_index| {
        series_array[sr_index] = SerieContainer.init(allocator);
        series_array[sr_index].serie = data[(sr_index * data_len)..((sr_index * data_len) + data_len)];
        series_array[sr_index].det_row = &det_mat[sr_index];
    }

    // var main_operator = DCCAOperator.init(time_steps[0..data_len], &series_array);

    var main_operator = allocator.create(DCCAOperator) catch unreachable;

    main_operator.time = time_steps[0..data_len];
    main_operator.previous = allocator.create(TimeData) catch unreachable;
    main_operator.previous.* = .{ .sum_x = 0, .sum_x2 = 0, .window_len = 0 };
    print("previous {}\n", .{main_operator.previous.window_len});

    main_operator.current = allocator.create(TimeData) catch unreachable;
    main_operator.current.* = .{ .sum_x = 0, .sum_x2 = 0, .window_len = 0 };
    print("current {}\n", .{main_operator.current.window_len});

    main_operator.series = series_array;

    // Time scale loop -- start
    for (0..tws_len, tws[0..]) |n_index, n| {
        defer {
            main_operator.*.reset();
            print("reset", .{});
        }
        const start_t_loop = Now_ms();
        print("\n---> starting calculations for time scale no. {} of size {}\n", .{
            n_index,
            n,
        });
        const start_alloc_time = Now_ms();
        // Allocation block -- start

        for (0..det_mat.len) |i| {
            det_mat[i] = allocator.alloc(f64, (n + 1)) catch |err| {
                print("Error allocating det_mat axis 1 line{}\n{}\n", .{ i, err });
                break;
            };
        }
        // Allocation block -- end
        const elapsed_alloc_time = Now_ms() - start_alloc_time;
        print("elapsed allocation time {}\n", .{elapsed_alloc_time});

        const no_of_windows: usize = (data_len - n);

        // for each Time window
        // time windows loop -- start

        const start_det_time = Now_ms();
        for (0..data_len - n) |win_start| {
            // print("wind_size: {} - window start: {}\n", .{ n, win_start });

            main_operator.shiftWindow(n, win_start, &F_DFA_arr[n_index]);

            for (0..data_count - 1) |pair_first| {
                for (pair_first + 1..data_count) |pair_sec| {
                    dcca_from_columns(det_mat, pair_first, pair_sec, n_index, &DCCA_arr);
                }
            }
        } // time windows loop -- end
        const elapsed_det_time = Now_ms() - start_det_time;
        print("elapsed det time {}\n", .{elapsed_det_time});

        //
        const start_filldfa_time = Now_ms();
        for (0..data_count) |sr_index| {
            F_DFA_arr[n_index][sr_index] = @sqrt(F_DFA_arr[n_index][sr_index] / @as(f64, @floatFromInt(no_of_windows)));
        }

        const elapsed_filldfa_time = Now_ms() - start_filldfa_time;
        print("elapsed fill dfa time {}\n", .{elapsed_filldfa_time});

        const start_dcca_pdcca_time = Now_ms();
        for (0..data_count - 1) |pair_first| {
            for (pair_first + 1..data_count) |pair_sec| {
                DCCA_arr[n_index][pair_first * data_count + pair_sec - 1] = DCCA_arr[n_index][
                    pair_first * data_count + pair_sec -
                        1
                ] / @as(f64, @floatFromInt(no_of_windows));

                if (P_dcca_matrix_output_flag == false) {
                    pdccaSimpleOutput(n_index, F_DFA_arr, DCCA_arr, pair_first, pair_sec, data_count,
                        // output
                        P_DCCA_arr);
                } else {
                    pdccaMatrixOutput(n_index, tws_len, data_count, F_DFA_arr, DCCA_arr, pair_first, pair_sec,
                        // output
                        P_DCCA_arr);
                }
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

export fn dcca(
    data: [*c]f64,
    data_len: usize,
    data_count: usize,
    tws: [*c]usize,
    tws_len: usize,
    time_steps: [*c]f64,
    //outputs
    F_DFA_arr: [*c][*c]f64,
    DCCA_arr: [*c][*c]f64,
) void {
    const start_time = Now_ms();
    print("Starting calculations:\n\nNumber of time steps:\t\t{}\nNumber of series:\t\t{}\nSize of series:\t\t\t{}\nNumber of DCCA combinations:\t{}\n\n", .{
        // initial values
        tws_len, data_count, data_len,
    });

    // Allocation setup -- start
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var arena = std.heap.ArenaAllocator.init(gpa.allocator());
    const allocator = arena.allocator();
    defer arena.deinit();
    // Allocation setup -- end

    var det_mat = allocator.alloc([]f64, data_count) catch unreachable;

    var series_array = allocator.alloc(*SerieContainer, data_count) catch unreachable;

    for (0..data_count) |sr_index| {
        series_array[sr_index] = SerieContainer.init(allocator);
        series_array[sr_index].serie = data[(sr_index * data_len)..((sr_index * data_len) + data_len)];
        series_array[sr_index].det_row = &det_mat[sr_index];
    }

    var main_operator = allocator.create(DCCAOperator) catch unreachable;

    main_operator.time = time_steps[0..data_len];
    main_operator.previous = allocator.create(TimeData) catch unreachable;
    main_operator.previous.* = .{ .sum_x = 0, .sum_x2 = 0, .window_len = 0 };
    print("previous {}\n", .{main_operator.previous.window_len});

    main_operator.current = allocator.create(TimeData) catch unreachable;
    main_operator.current.* = .{ .sum_x = 0, .sum_x2 = 0, .window_len = 0 };
    print("current {}\n", .{main_operator.current.window_len});

    main_operator.series = series_array;

    // Time scale loop -- start
    for (0..tws_len, tws[0..]) |n_index, n| {
        defer {
            main_operator.*.reset();
            print("reset", .{});
        }
        const start_t_loop = Now_ms();
        print("\n---> starting calculations for time scale no. {} of size {}\n", .{
            n_index,
            n,
        });
        const start_alloc_time = Now_ms();
        // Allocation block -- start

        for (0..det_mat.len) |i| {
            det_mat[i] = allocator.alloc(f64, (n + 1)) catch |err| {
                print("Error allocating det_mat axis 1 line{}\n{}\n", .{ i, err });
                break;
            };
        }
        // Allocation block -- end
        const elapsed_alloc_time = Now_ms() - start_alloc_time;
        print("elapsed allocation time {}\n", .{elapsed_alloc_time});

        const no_of_windows: usize = (data_len - n);

        // for each Time window
        // time windows loop -- start

        const start_det_time = Now_ms();
        for (0..data_len - n) |win_start| {
            // print("wind_size: {} - window start: {}\n", .{ n, win_start });

            main_operator.shiftWindow(n, win_start, &F_DFA_arr[n_index]);

            for (0..data_count - 1) |pair_first| {
                for (pair_first + 1..data_count) |pair_sec| {
                    dcca_from_columns(det_mat, pair_first, pair_sec, n_index, &DCCA_arr);
                }
            }
        } // time windows loop -- end
        const elapsed_det_time = Now_ms() - start_det_time;
        print("elapsed det time {}\n", .{elapsed_det_time});

        //
        const start_filldfa_time = Now_ms();
        for (0..data_count) |sr_index| {
            F_DFA_arr[n_index][sr_index] = @sqrt(F_DFA_arr[n_index][sr_index] / @as(f64, @floatFromInt(no_of_windows)));
        }

        const elapsed_filldfa_time = Now_ms() - start_filldfa_time;
        print("elapsed fill dfa time {}\n", .{elapsed_filldfa_time});

        const start_dcca_pdcca_time = Now_ms();
        for (0..data_count - 1) |pair_first| {
            for (pair_first + 1..data_count) |pair_sec| {
                DCCA_arr[n_index][pair_first * data_count + pair_sec - 1] = DCCA_arr[n_index][
                    pair_first * data_count + pair_sec - 1
                ] / @as(f64, @floatFromInt(no_of_windows));
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

pub fn dcca_from_columns(
    det_mat: [][]f64,
    DCCA_of: *const [*c][*c]usize,
    dcca_pair: usize,
    n_index: usize,

    // output
    DCCA_arr: *const [*c][*c]f64,
) void {
    var dcca_temp: f64 = 0;
    for (0..det_mat[0].len) |j| {
        dcca_temp += det_mat[DCCA_of.*[dcca_pair][0]][j] * det_mat[DCCA_of.*[dcca_pair][1]][j];
    }
    DCCA_arr.*[n_index][dcca_pair] += dcca_temp / @as(f64, @floatFromInt(det_mat[0].len));
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
    P_DCCA_arr[n_index * DCCA_of_n_rows + dcca_pair] = DCCA_arr[n_index][dcca_pair] /
        (F_DFA_arr[n_index][DCCA_of[dcca_pair][0]] * F_DFA_arr[n_index][DCCA_of[dcca_pair][1]]);
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
