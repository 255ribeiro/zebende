const std = @import("std");
const pow = std.math.pow;
const Thread = std.Thread;
const Mutex = std.Thread.Mutex;
const print = std.debug.print;
const expect = std.testing.expect;

// c import

// const c = @cImport({
//     // @cInclude("stdio.h");
//     @cInclude("stdlib.h");
// });

// const testing = std.testing;

const DfaDccaOperator = struct {
    // pre seted
    lock: Mutex,
    count: usize,
    accumulator: f64,

    fn init(allocator: std.mem.Allocator) *DfaDccaOperator {
        const ptr = allocator.create(DfaDccaOperator) catch unreachable;
        ptr.count = 0;
        ptr.accumulator = 0;
        ptr.lock.unlock();
        return ptr;
    }

    fn calcDfa(self: *DfaDccaOperator, det_mat: *const [][][]f64, win_index: usize, sr_index: usize, n: usize, no_of_windows: usize, F_DFA_arr_ptr: [*c]f64) void {
        var temp: f64 = 0;
        for (0..(n + 1)) |i| {
            temp += pow(f64, det_mat.*[i][win_index][sr_index], 2);
        }

        // lock accumulator and counter
        self.lock.lock();
        defer self.lock.unlock();
        // write values
        self.accumulator += temp / @as(f64, @floatFromInt(n + 1));
        self.count += 1;
        // stop condition
        if (self.count == no_of_windows) {
            F_DFA_arr_ptr[sr_index] += @sqrt(self.accumulator / @as(f64, @floatFromInt(no_of_windows)));
            print("DFA calculated for serie no. {d} for {d} time widwos\n", .{ sr_index, self.count });
        }
    }

    // DCCA calc function
    fn calcDcca(
        self: *DfaDccaOperator,
        det_mat: *const [][][]f64,
        win_index: usize,
        DCCA_of_ptr: *const [*c][*c]usize,
        dcca_pair: usize,
        n: usize,
        no_of_windows: usize,
        // output
        DCCA_arr_ptr: [*c]f64,
    ) void {
        var temp: f64 = 0;
        for (0..(n + 1)) |i| {
            temp += det_mat.*[i][win_index][DCCA_of_ptr.*[dcca_pair][0]] * det_mat.*[i][win_index][DCCA_of_ptr.*[dcca_pair][1]];
        }

        // lock accumulator and counter
        self.lock.lock();
        defer self.lock.unlock();
        // write values
        self.accumulator += temp / @as(f64, @floatFromInt(n + 1));
        self.count += 1;
        // stop condition
        if (self.count == no_of_windows) {
            DCCA_arr_ptr[dcca_pair] += self.accumulator / @as(f64, @floatFromInt(no_of_windows));
            print("DCCA calculated for DCCA of no. {d} for {d} windows\n", .{ dcca_pair, self.count });
        }
    }
};

pub fn srData(data_index: usize, data_count: usize) usize {
    const sr_index = data_index % data_count;
    return sr_index;
}

pub fn winData(data_index: usize, data_count: usize) usize {
    const win_index = @divFloor(data_index, data_count);
    return win_index;
}

// Dfa warpper function
pub fn dfaWriter(
    dfa_operator: *DfaDccaOperator,
    det_mat: *const [][][]f64,
    win_index: usize,
    sr_index: usize,
    n: usize,
    no_of_windows: usize,
    // output
    F_DFA_arr_ptr: [*c]f64,
) void {
    dfa_operator.*.calcDfa(det_mat, win_index, sr_index, n, no_of_windows,
    // output
    F_DFA_arr_ptr);
}

pub fn dccaWriter(
    dcca_operator: *DfaDccaOperator,
    det_mat: *const [][][]f64,
    win_index: usize,
    DCCA_of_ptr: *const [*c][*c]usize,
    dcca_pair: usize,
    n: usize,
    no_of_windows: usize,
    // output
    DCCA_arr_ptr: [*c]f64,
) void {
    dcca_operator.*.calcDcca(det_mat, win_index, DCCA_of_ptr, dcca_pair, n, no_of_windows,
    // output
    DCCA_arr_ptr);
}

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
    const no_of_threads = Thread.getCpuCount() catch unreachable;
    print("Starting calculations:\n\nNumber of time steps:\t\t{}\nNumber of series:\t\t{}\nSize of series:\t\t\t{}\nNumber of DCCA combinations:\t{}\nNumber of threads:\t\t{}\n", .{
        // initial values
        tws_len, data_count, data_len, DCCA_of_n_rows, no_of_threads,
    });

    // Allocation setup -- start
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    // Allocation setup -- end

    // Time scale loop -- start
    for (0..tws_len) |n_index| {
        const n = tws[n_index];
        const no_of_windows: usize = (data_len - n);
        print("\n---> starting calculations for time scale no. {} of size {} - number of windows: {}\n", .{ n_index, n, no_of_windows });

        // Allocation block
        print("Allocating {} f64 in detrended matrix\n", .{(n + 1) * no_of_windows * data_count});
        var arena = std.heap.ArenaAllocator.init(gpa.allocator());
        const allocator = arena.allocator();
        defer arena.deinit();

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

        { // handlers detrend scope --start
            var hand_det = allocator.alloc(std.Thread, (no_of_windows * data_count)) catch |err| {
                print("Error allocating handels {}\n", .{err});
                break;
            };

            defer allocator.free(hand_det);

            // windows loop -- start
            for (0..(no_of_windows * data_count)) |data_index| {

                // spawn hadlers
                hand_det[data_index] = std.Thread.spawn(.{}, detendPrl, .{
                    &data,
                    &time_steps,
                    // data index
                    data_index,
                    // time scale args
                    n,
                    no_of_windows,
                    // mat output
                    &det_mat,
                }) catch |err| {
                    print("Error spawning thread hadels no. {} error: {}\n", .{ data_index, err });
                    break;
                };
            } // windows loop -- end

            print("detrend handlers in use {d} from a total of {d} at time scale no. {d} of size {d}\n", .{ hand_det.len, (no_of_windows * data_count), n_index, n });

            for (hand_det) |h| h.join();

            print("detrend threads joined for time scale no. {} of size {} \n", .{ n_index, n });
        } // handlers detrend scope -- end

        // Handlers DFA DCCA scope -- start

        {
            var operator_alloc = allocator.alloc(*DfaDccaOperator, (data_count + DCCA_of_n_rows)) catch |err| {
                print("Error allocating dfa/dcca operators\n{}\n", .{err});
                break;
            };

            defer allocator.free(operator_alloc);

            for (0..operator_alloc.len) |i| {
                operator_alloc[i] = DfaDccaOperator.init(allocator);
                operator_alloc[i].lock.unlock();
            }

            var hand_dfa_dcca = allocator.alloc(std.Thread, no_of_windows * operator_alloc.len) catch |err| {
                print("Error allocating handels for dfa/dcca operators {}\n", .{err});
                break;
            };

            defer allocator.free(hand_dfa_dcca);

            print("DFA/DCCA operators initialized successfully for time scale no. {} of size {}\n", .{ n_index, n });

            for (0..hand_dfa_dcca.len) |data_index| {
                const op_index = srData(data_index, (data_count + DCCA_of_n_rows));
                const win_index = winData(data_index, (data_count + DCCA_of_n_rows));

                if (op_index < data_count) { // dfa calculations
                    hand_dfa_dcca[data_index] = std.Thread.spawn(.{}, dfaWriter, .{
                        operator_alloc[op_index],
                        // read only matrix
                        &det_mat,
                        // window id
                        win_index,
                        //series id
                        op_index,
                        // time scale args
                        n,
                        no_of_windows,
                        // final output
                        &F_DFA_arr[n_index].*,
                    }) catch |err| {
                        print("err: {}", .{err});
                        break;
                    };
                } else {
                    hand_dfa_dcca[data_index] = std.Thread.spawn(.{}, dccaWriter, .{
                        operator_alloc[op_index],
                        //read only matrix
                        &det_mat,
                        // win index
                        win_index,
                        // DCCA args
                        &DCCA_of,
                        (op_index - data_count),
                        // time scale args
                        n,
                        no_of_windows,
                        // output
                        &DCCA_arr[n_index].*,
                    }) catch |err| {
                        print("err: {}", .{err});
                        break;
                    };
                }
            }

            for (hand_dfa_dcca) |h| h.join();

            print("DFA and DCCA calculated for time scale no. {} of size {}\n", .{ n_index, n });
        } // handles DFA DCCA scope -- end

        // // p_dcca calc

        var hand_pdcca = allocator.alloc(std.Thread, DCCA_of_n_rows) catch unreachable;
        for (0..DCCA_of_n_rows) |dcca_pair| {
            hand_pdcca[dcca_pair] = Thread.spawn(.{}, pDccaSimpleOutPrl, .{
                // mat imputs
                &F_DFA_arr[n_index].*,
                &DCCA_arr[n_index].*,
                &DCCA_of,
                // usize paris
                dcca_pair,
                // output
                &P_DCCA_arr[n_index].*,
            }) catch unreachable;
        }

        for (hand_pdcca) |h| h.join();

        print("P_DCCA calculated for time scale no. {} of size {}\n", .{ n_index, n });

        // Allocation free
    } // Time scale loop -- end
} // function --end

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

pub fn detendPrl(
    // data input
    data: *const [*c]f64,
    time_steps: *const [*c]f64,
    // data index
    data_index: usize,
    // time scale args
    n: usize,
    no_of_windows: usize,
    // mat output
    det_mat: *const [][][]f64,
) void {
    const tw_start: usize = data_index % no_of_windows;
    const sr_ind = @divFloor(data_index, no_of_windows);
    const sw_start: usize = data_index + (sr_ind * n);

    const time = time_steps.*[tw_start .. tw_start + (n + 1)];
    const win = data.*[sw_start .. sw_start + (n + 1)];
    const slp_int = linLestSquaresFit(win, time);
    for (win, time, 0..) |w, t, i| {
        det_mat.*[i][tw_start][sr_ind] = w - (slp_int[0] * t + slp_int[1]);
    }
}

pub fn pDccaSimpleOutPrl(
    // mat imputs
    F_DFA_arr_ptr: [*c]f64,
    DCCA_arr_ptr: [*c]f64,
    // dcca of
    DCCA_of: *const [*c][*c]usize,
    // usize vars
    dcca_pair: usize,
    // output
    P_DCCA_arr_ptr: [*c]f64,
) void {
    P_DCCA_arr_ptr[dcca_pair] = DCCA_arr_ptr[dcca_pair] / (F_DFA_arr_ptr[DCCA_of.*[dcca_pair][0]] * F_DFA_arr_ptr[DCCA_of.*[dcca_pair][1]]);
}
