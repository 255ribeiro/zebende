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

const SeriesData = struct {
    sum_x: f64 = 0,
    sum_x2: f64 = 0,
    window_len: usize = 0,
    sum_y: []f64,
    sum_xy: []f64,

    fn init(allocator: std.mem.Allocator, data_count: usize) *SeriesData {
        const ptr = allocator.create(SeriesData) catch @panic("error allocating series data");
        ptr.sum_x = 0;
        ptr.sum_x2 = 0;
        ptr.window_len = 0;

        ptr.sum_y = allocator.alloc(f64, data_count) catch @panic("error allocating sum y");

        ptr.sum_xy = allocator.alloc(f64, data_count) catch @panic("error allocating sum xy");

        for (0..data_count) |i| {
            ptr.sum_y[i] = 0;
            ptr.sum_xy[i] = 0;
        }

        return ptr;
    }
};

const MainOperator = struct {
    left_x: f64,
    previous: *SeriesData,
    current: *SeriesData,
    left_y: []f64,
    slope: []f64,
    inter: []f64,
    dfa_temp: []f64,
    time: []f64,
    time_window: []f64,
    window_pointer: [*]f64,
    det_mat: [][]f64,
    data: [*]f64,
    allocator: std.mem.Allocator,

    fn init(allocator: std.mem.Allocator, data_count: usize, time: []f64) *MainOperator {
        const ptr = allocator.create(MainOperator) catch unreachable;
        ptr.left_x = 0;
        ptr.allocator = allocator;
        ptr.previous = SeriesData.init(allocator, data_count);
        ptr.current = SeriesData.init(allocator, data_count);

        ptr.left_y = allocator.alloc(f64, data_count) catch unreachable;

        ptr.slope = allocator.alloc(f64, data_count) catch unreachable;

        ptr.inter = allocator.alloc(f64, data_count) catch unreachable;

        ptr.dfa_temp = allocator.alloc(f64, data_count) catch unreachable;

        ptr.time = time;

        for (0..data_count) |i| {
            ptr.left_y[i] = 0;
            ptr.slope[i] = 0;
            ptr.inter[i] = 0;
            ptr.dfa_temp[i] = 0;
        }

        return ptr;
    }
    // shift window
    fn shiftWindow(self: *MainOperator, n: usize, win_start: usize, F_DFA_ptr: [*c]f64) void {
        self.time_window = self.time[win_start..][0..(n + 1)];

        if (win_start != 0) {
            self.current.sum_x = self.current.sum_x - self.left_x + self.time_window[n];
            self.current.sum_x2 = self.current.sum_x2 - pow(f64, self.left_x, 2) + pow(f64, self.time_window[n], 2);
            // series window pointer -- set
            self.window_pointer = self.data + ((win_start + n) * self.left_y.len);

            for (0..self.left_y.len) |sr_index| {
                // series window pointer -- increment
                self.window_pointer += sr_index;
                self.current.sum_y[sr_index] = self.current.sum_y[sr_index] - self.left_y[sr_index] + self.window_pointer[0];
                self.current.sum_xy[sr_index] = self.current.sum_xy[sr_index] - (self.left_x * self.left_y[sr_index]) + self.window_pointer[0];
            }
        } else { // ( win_start == 0)

            print("else \n", .{});

            self.current.window_len = n + 1;

            for (self.previous.window_len..self.current.window_len) |i| {
                self.previous.sum_x += self.time_window[i];

                self.previous.sum_x2 += pow(f64, self.time_window[i], 2);
                // series window pointer -- set
                self.window_pointer = self.data + (i * self.left_y.len);

                for (0..self.left_y.len) |sr_index| {
                    //series window pointer -- increment
                    self.window_pointer += sr_index;
                    self.previous.sum_y[sr_index] += self.window_pointer[0];
                    self.previous.sum_xy[sr_index] += self.time_window[i] * self.window_pointer[0];
                }
            }

            // updating current sum values
            // x and x^2 sums
            self.current.sum_x = self.previous.sum_x;
            self.current.sum_x2 = self.previous.sum_x2;
            // y and xy sums
            @memcpy(self.current.sum_y, self.previous.sum_y);
            @memcpy(self.current.sum_xy, self.previous.sum_xy);
        }
        // slope inter calculations
        const len: f64 = @as(f64, @floatFromInt(self.time_window.len));

        for (0..self.left_y.len) |sr_index| {
            //slope
            self.slope[sr_index] = (((len * self.current.sum_xy[sr_index]) - (self.current.sum_x * self.current.sum_y[sr_index])) /
                ((len * self.current.sum_x2) - (pow(f64, self.current.sum_x, 2))));
            //inter
            self.inter[sr_index] = ((self.current.sum_y[sr_index] - (self.slope[sr_index] * self.current.sum_x)) /
                (len));

            // series window pointer -- set
            // self.window_pointer = self.data + (win_start * self.left_y.len) + sr_index;
            // reset dfa_temp
            self.dfa_temp[sr_index] = 0;
            for (self.time_window, 0..) |t, j| {
                //series window pointer -- increment
                self.window_pointer = self.data + ((win_start + j) * self.left_y.len);

                self.det_mat[sr_index][j] = self.window_pointer[sr_index] - (self.slope[sr_index] * t + self.inter[sr_index]);
                self.dfa_temp[sr_index] += pow(f64, self.det_mat[sr_index][j], 2);
            }
            F_DFA_ptr[sr_index] += self.dfa_temp[sr_index] / len;
        }

        // detrended and dfa calculations

        // reset left values
        self.left_x = self.time_window[0];
        self.left_y = self.data[win_start * self.left_y.len ..][0..self.left_y.len];
    }

    // reset function
    fn reset(
        self: MainOperator,
    ) void {
        self.previous.window_len = self.current.window_len;
    }
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
    var arena = std.heap.ArenaAllocator.init(gpa.allocator());
    const allocator = arena.allocator();
    defer arena.deinit();
    // Allocation setup -- end

    var det_mat = allocator.alloc([]f64, data_count) catch unreachable;

    // var main_operator = MainOperator.init(time_steps[0..data_len], &series_array);

    var main_operator = allocator.create(MainOperator) catch unreachable;

    main_operator = MainOperator.init(allocator, data_count, time_steps[0..data_len]);

    main_operator.data = data;
    main_operator.det_mat = det_mat;

    main_operator.previous = SeriesData.init(allocator, data_count);
    // main_operator.previous.* = .{ .sum_x = 0, .sum_x2 = 0, .window_len = 0, .sum_y = .{0} ** data_count };
    print("previous {d}\n", .{main_operator.left_y});

    main_operator.current = SeriesData.init(allocator, data_count);
    // main_operator.current.* = .{ .sum_x = 0, .sum_x2 = 0, .window_len = 0, .sum_y = .{0} ** data_count };
    print("current {d}\n", .{main_operator.current.window_len});

    // Time scale loop -- sart
    for (0..tws_len, tws[0..]) |n_index, n| {
        defer {
            main_operator.reset();
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

            main_operator.shiftWindow(n, win_start, F_DFA_arr[n_index]);
            // print("shifted!\n", .{});

            for (0..DCCA_of_n_rows) |dcca_pair| {
                dcca(det_mat, DCCA_of, dcca_pair, n_index, DCCA_arr);
            }

            // print("time window end\n", .{});
        } // time windows loop -- end
        const elapsed_det_time = Now_ms() - start_det_time;
        print("elapsed det time {}\n", .{elapsed_det_time});

        //
        const start_filldfa_time = Now_ms();
        print("fill dfa \n", .{});

        for (0..data_count) |sr_index| {
            F_DFA_arr[n_index][sr_index] = @sqrt(F_DFA_arr[n_index][sr_index] / @as(f64, @floatFromInt(no_of_windows)));
            print("F_DFA_arr[n_index][sr_index] {d}\n", .{F_DFA_arr[n_index][sr_index]});
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
        // Allocation free
        const tws_elapsed = (Now_ms() - start_t_loop);
        print("tws elapsed: {} milis for time scale no. {} of size {}\n", .{ tws_elapsed, n_index, n });
    } // time window scale -- loop end

    print("elapsed time for all {}\n", .{Now_ms() - start_time});
} // function -- end

pub fn dcca(
    det_mat: [][]f64,
    DCCA_of: [*c][*c]usize,
    dcca_pair: usize,
    n_index: usize,
    // output
    DCCA_arr: [*c][*c]f64,
) void {
    var dcca_temp: f64 = 0;
    for (0..det_mat[0].len) |j| {
        dcca_temp += det_mat[DCCA_of[dcca_pair][0]][j] * det_mat[DCCA_of[dcca_pair][1]][j];
    }

    DCCA_arr[n_index][dcca_pair] += dcca_temp / @as(f64, @floatFromInt(det_mat[0].len));
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
