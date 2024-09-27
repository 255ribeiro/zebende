import numpy as np

def arr_2d_to_c(arr):
    return (arr.__array_interface__['data'][0] + np.arange(arr.shape[0] ) * arr.strides[0]).astype(np.uintp)