import numpy as np
from typing import Any
from numpy.typing import NDArray

def arr_2d_to_c(arr:NDArray[Any])-> NDArray[Any]:
    """Convert a 2D numpy array to a array of pointers that can be passed to a C ABI funtion.
    Args:
        arr (NDArray[Any]): 2D array of any type.
    Returns:
        NDArray[Any]: 1D array of pointers.
    """
        # Ensure the array is C-contiguous
    arr = np.ascontiguousarray(arr)
    
    # Safely access the base address
    base_address = arr.__array_interface__['data']
    if isinstance(base_address, tuple):  # Handle older NumPy versions
        base_address = base_address[0]
    
    # Compute row pointers
    row_pointers = (base_address + np.arange(arr.shape[0]) * arr.strides[0]).astype(np.uintp)
    
    return np.ascontiguousarray(row_pointers)