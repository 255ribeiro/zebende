import numpy as np


def arr_2d_to_c(arr):
        # Ensure the array is C-contiguous
    arr = np.ascontiguousarray(arr)
    
    # Safely access the base address
    base_address = arr.__array_interface__['data']
    if isinstance(base_address, tuple):  # Handle older NumPy versions
        base_address = base_address[0]
    
    # Compute row pointers
    row_pointers = (base_address + np.arange(arr.shape[0]) * arr.strides[0]).astype(np.uintp)
    
    return np.ascontiguousarray(row_pointers)