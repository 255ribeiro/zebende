const std = @import("std");
const print = std.debug.print;
const expect = std.testing.expect;
const pow = std.math.pow;
const Now_ms = std.time.microTimestamp;
// const testing = std.testing;

// c import
// const c = @cImport({
//     // @cInclude("stdio.h");
//     @cInclude("stdlib.h");
// });

const Seriedata = struct {
    sum_x: f64 = 0,
    sum_x2: f64 = 0,
    sum_y: f64 = 0,
    sum_xy: f64 = 0,

    window_size: usize = 0,
    fn init(allocator: std.mem.Allocator) *Seriedata {
        const ptr = allocator.create(Seriedata) catch unreachable;
        ptr.sum_x = 0;
        ptr.sum_x2 = 0;
        ptr.sum_xy = 0;
        ptr.sum_y = 0;
        ptr.window_size = 0;
        return ptr;
    }
};

const SerieContainer = struct {
    previous: *Seriedata,
    current: *Seriedata,
    serie: []f64,
    det_row: *const []f64,
    left_x: f64,
    left_y: f64,
    slope: f64,
    inter: f64,
    dfa_temp: f64,

    fn init(allocator: std.mem.Allocator) *SerieContainer {
        const ptr = allocator.create(SerieContainer) catch unreachable;
        ptr.previous = Seriedata.init(allocator);
        ptr.current = Seriedata.init(allocator);
        return ptr;
    }
    fn shiftWindow(self: *SerieContainer, n: usize, win_start: usize, time: []f64, F_DFA_ptr: *allowzero f64) void {
        const window = self.serie[win_start .. win_start + (n + 1)];
        if (self.previous.window_size != self.current.window_size) {
            // print("defined win start {}\n", .{win_start});
            self.current.sum_x = self.current.sum_x - self.left_x + time[time.len - 1];
            self.current.sum_x2 = self.current.sum_x2 - pow(f64, self.left_x, 2) + pow(f64, time[time.len - 1], 2);
            self.current.sum_y = self.current.sum_y - self.left_y + window[window.len - 1];
            self.current.sum_xy = self.current.sum_xy - (self.left_x * self.left_y) + (time[time.len - 1] * window[window.len - 1]);
        } else {
            self.current.window_size = n + 1;
            // print("undefined win start {}: {d:0.10} time:{d}\n", .{ win_start, self.previous.sum_x, time });
            for (self.previous.window_size..self.current.window_size) |i| {
                self.previous.sum_x += time[i];
                self.previous.sum_x2 += pow(f64, time[i], 2);
                self.previous.sum_y += window[i];
                self.previous.sum_xy += time[i] * window[i];
            }
            self.current.sum_x = self.previous.sum_x;
            self.current.sum_x2 = self.previous.sum_x2;
            self.current.sum_xy = self.previous.sum_xy;
            self.current.sum_y = self.previous.sum_y;
        }
        self.left_x = time[0];
        self.left_y = window[0];

        const len: f64 = @as(f64, @floatFromInt(self.current.window_size));
        //slope
        self.slope = (((len * self.current.sum_xy) - (self.current.sum_x * self.current.sum_y)) /
            ((len * self.current.sum_x2) - (pow(f64, self.current.sum_x, 2))));
        //inter
        self.inter = ((self.current.sum_y - (self.slope * self.current.sum_x)) /
            (len));
        //result
        self.dfa_temp = 0;
        for (window, time, 0..) |w, t, i| {
            self.det_row.*[i] = w - (self.slope * t + self.inter);
            self.dfa_temp += pow(f64, self.det_row.*[i], 2);
        }
        F_DFA_ptr.* += self.dfa_temp / len;
    }

    fn reset(
        self: *SerieContainer,
    ) void {
        self.previous.window_size = self.current.window_size;
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

    var series_array = allocator.alloc(*SerieContainer, data_count) catch unreachable;

    for (0..data_count) |sr_index| {
        series_array[sr_index] = SerieContainer.init(allocator);
        series_array[sr_index].serie = data[(sr_index * data_len)..((sr_index * data_len) + data_len)];
        series_array[sr_index].det_row = &det_mat[sr_index];

        print("ptr {d} : sum x previous {d} sum x current {d}\n", .{ sr_index, series_array[sr_index].previous.sum_x, series_array[sr_index].current.sum_x });
    }

    // Time scale loop -- sart
    for (0..tws_len, tws[0..]) |n_index, n| {
        defer {
            for (series_array) |serie| serie.reset();
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

            // time windows slice
            const time_window = time_steps[win_start .. win_start + (n + 1)];

            // for each data series
            // series loop -- start
            for (series_array, 0..) |serie, sr_index| {
                serie.shiftWindow(n, win_start, time_window, &F_DFA_arr[n_index][sr_index]);
            } // series loop -- end

            for (0..DCCA_of_n_rows) |dcca_pair| {
                dcca_from_columns(det_mat, &DCCA_of, dcca_pair, n_index, &DCCA_arr);
            }
        } // time windows loop -- end
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