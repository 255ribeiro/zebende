import numpy as np


# linaer Least Squares fit
def liner_ls_fit(arr_x, arr_y):
    n = arr_x.size

    x_sum = 0
    y_sum = 0
    xy_sum = 0
    x2_sum = 0
    for i in range(n):
        x_sum += arr_x[i]
        y_sum += arr_y[i,:]
        xy_sum += arr_x[i] * arr_y[i,:]
        x2_sum += arr_x[i]**2
    # slope (m in y = mx +b)
    slope = (
            ((n * xy_sum) - (x_sum * y_sum)) /
            ((n * x2_sum) - (x_sum**2))
            )
    # intercept (b in y = mx +b)
    inter = (
            (y_sum - (slope * x_sum)) /
            (n)
            )

    return slope, inter


# detrended series
def detrended_series(arr_x, arr_y):
    n = arr_x.size
    slope, inter = liner_ls_fit(arr_x, arr_y)
    out = np.empty_like(arr_y)
    for i in range(n):
        out[i] = arr_y[i,:] - (slope * arr_x[i] + inter)
    return out
