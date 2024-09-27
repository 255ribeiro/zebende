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

const Dfa_struct = struct {
    // pre seted
    lock: Mutex,
    count: usize,
    dfa_accumulator: f64,

    fn init(allocator: std.mem.Allocator) *Dfa_struct {
        const ptr = allocator.create(Dfa_struct) catch unreachable;
        ptr.count = 0;
        ptr.dfa_accumulator = 0;
        ptr.lock.unlock();
        return ptr;
    }

    fn calcDfa(self: *Dfa_struct, det_mat: *const [][][]f64, win_index: usize, sr_index: usize, n: usize, no_of_windows: usize, F_DFA_arr_ptr: [*c]f64) void {
        var dfa_temp: f64 = 0;
        for (0..(n + 1)) |i| {
            dfa_temp += pow(f64, det_mat.*[i][win_index][sr_index], 2);
        }

        // lock accumulator and counter
        self.lock.lock();
        defer self.lock.unlock();
        // write values
        self.dfa_accumulator += dfa_temp / @as(f64, @floatFromInt(n + 1));
        self.count += 1;
        // stop condition
        if (self.count == no_of_windows) {
            F_DFA_arr_ptr[sr_index] += @sqrt(self.dfa_accumulator / @as(f64, @floatFromInt(no_of_windows)));
            print("dfa series {d} count {d}\n", .{ sr_index, self.count });
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

pub fn dfaWriter(dfa_struct: *Dfa_struct, det_mat: *const [][][]f64, win_index: usize, sr_index: usize, n: usize, no_of_windows: usize, F_DFA_arr_ptr: [*c]f64) void {
    dfa_struct.*.calcDfa(det_mat, win_index, sr_index, n, no_of_windows,
    // output
    F_DFA_arr_ptr);
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
    _ = DCCA_of;
    _ = DCCA_of_n_rows;

    _ = DCCA_arr;

    print(" tipe of {}\n", .{@TypeOf(P_DCCA_arr)});

    // Allocation setup -- start
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    // Allocation setup -- end

    // Time scale loop -- start
    for (0..tws_len) |n_index| {
        const n = tws[n_index];
        const no_of_windows: usize = (data_len - n);

        // Allocation block

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

            // windows loop -- start
            for (0..(no_of_windows * data_count)) |data_index| {

                // spawn hadlers
                var det_threads = std.Thread.spawn(.{}, detend_prl, .{
                    &data,
                    &time_steps,
                    // usize vars
                    data_index,
                    no_of_windows,
                    n,
                    // mat output
                    &det_mat,
                }) catch |err| {
                    print("Error spawning thread hadels no. {} error: {}\n", .{ data_index, err });
                    break;
                };

                det_threads.join();
            } // windows loop -- end

        } // handlers detrend scope -- end

        // // dfa calc

        {
            var dfa_alloc = allocator.alloc(*Dfa_struct, data_count) catch |err| {
                print("Error allocating dfa array\n{}\n", .{err});
                break;
            };

            for (0..dfa_alloc.len) |i| {
                dfa_alloc[i] = Dfa_struct.init(allocator);
                dfa_alloc[i].lock.unlock();

                print("dfa count: {d}\ndfa accumulator: {d} {}\n", .{ dfa_alloc[i].count, dfa_alloc[i].dfa_accumulator, dfa_alloc[i].lock.impl });
            }

            for (0..(no_of_windows * data_count)) |data_index| {
                var thread_dfa = std.Thread.spawn(.{}, dfaWriter, .{
                    dfa_alloc[srData(data_index, data_count)],
                    // reading matrix
                    &det_mat,
                    // window id
                    winData(data_index, data_count),
                    //series id
                    srData(data_index, data_count),
                    // args
                    n,
                    no_of_windows,
                    // final output
                    &F_DFA_arr[n_index].*,
                }) catch |err| {
                    print("err: {}", .{err});
                    break;
                };
                thread_dfa.join();
            }
        } // handles F_DFA scope -- end

        // // dcca calc

        // { // poll DCCA scope -- start
        //     var pool: std.Thread.Pool = undefined;
        //     pool.init(.{ .allocator = allocator }) catch |err| {
        //         print("Error allocating thread pool {}\n", .{err});
        //         break;
        //     };
        //     defer pool.deinit();

        //     for (0..(no_of_windows * DCCA_of_n_rows)) |data_index| {
        //         pool.spawn(DCCA_fill_prl, .{
        //             // input matrices
        //             &det_mat,   &DCCA_of,
        //             // usize vars
        //             data_index, n,
        //             n_index,    DCCA_of_n_rows,
        //             // aoutput mat
        //             &DCCA_arr,
        //         }) catch |err| {
        //             print("Error spawning thread pool no. {} error: {}\n", .{ data_index, err });
        //             break;
        //         };
        //     }
        // } // poll DCCA scope -- end

        // for (0..DCCA_of_n_rows) |dcca_pair| {
        //     DCCA_arr[n_index][dcca_pair] = DCCA_arr[n_index][dcca_pair] / @as(f64, @floatFromInt(no_of_windows));
        // }
        // // p_dcca calc

        // for (0..DCCA_of_n_rows) |dcca_pair| {
        //     DCCA_arr[n_index][dcca_pair] = DCCA_arr[n_index][dcca_pair] / @as(f64, @floatFromInt(no_of_windows));

        //     p_dcca_simple_output(
        //     // mat imputs
        //     &F_DFA_arr, &DCCA_arr, &DCCA_of,
        //     // usize paris
        //     dcca_pair, n_index,
        //     // output
        //     &P_DCCA_arr);
        // }

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
