import ctypes

import numpy as np
from typing import Literal

from .get_platform import get_platform

from .array_to_c_pointer_convert import arr_2d_to_c

from numpy.typing import NDArray
from numpy.ctypeslib import ndpointer

uint_c_type = np.uintp

ENUM_DCCA_of = Literal['all']

def dfa(
        input_data: NDArray[np.float64], 
        tws:  NDArray[np.int64] | NDArray[np.float64],
        axis: int = 0          
        ) ->    tuple[
                NDArray[np.float64],    # DFA
                ]:

    """A function that calculates the 
        DFA (Detrended Fluctuation Analysis)
        for a group of time series

        Args:
            input_data (NDArray[np.float64]): 2D array of times series integrated data.
            tws (NDArray[np.float64] | NDArray[np.float64]): 1D array of time scales.
            axis (int): axis of the input_data matrix that contains the values for each time series. Defalts to 0.

        Returns:
            DFA(NDArray[np.float64]):_description_. 
           
    """
    
    assert (tws[:-1] < tws[1:]).all() == True , ("""time window scales (tws) values must be in crescent order.""")


    # setting lib path

    zz = ctypes.cdll.LoadLibrary(get_platform())

    # outputs types
    c_2d_any_1d_uint = ndpointer(dtype = uint_c_type, ndim=1, flags='C')
    c_1d_double = ndpointer(dtype=np.double, ndim=1, flags='C')
    # getting the module file
    _dfa = zz.dfa
    # setting outputs
    _dfa.argtypes = [ 
                    c_1d_double, ctypes.c_size_t, ctypes.c_size_t, # input data
                    c_2d_any_1d_uint, ctypes.c_size_t, # tws
                    c_1d_double, # Time steps
                    c_2d_any_1d_uint, # DFA_arr
                    ]

    # Preparing input data    

    if not tws.flags.c_contiguous:
        tws = np.ascontiguousarray(tws)

    
    # preparing tws array
    tws = tws.astype(uint_c_type)

    # preparing output array
    F_DFA_arr = np.ascontiguousarray(np.zeros(shape=(tws.shape[0], input_data.shape[axis]), dtype=input_data.dtype))
    c_DFA_arr = arr_2d_to_c(F_DFA_arr)

    # preparing data
    data_shape = input_data.shape

    if input_data.ndim == 1:   
        x_len = input_data.size
        x_cnt = 1
        input_data = np.ascontiguousarray(input_data.flatten())

    elif axis == 0:
        x_len = data_shape[1]
        x_cnt = data_shape[0]
        input_data = np.ascontiguousarray(input_data.flatten())

    elif axis == 1:
        x_len = data_shape[0]
        x_cnt = data_shape[1]
        input_data = np.ascontiguousarray(np.asfortranarray(input_data).flatten(order="K"))

    time_steps = np.ascontiguousarray(np.arange(x_len, dtype=input_data.dtype))

    # calling functon  
    _dfa(input_data,  x_len, x_cnt, 
            tws, tws.size, 
            time_steps, 
            c_DFA_arr
            )

    return F_DFA_arr
