import platform
import ctypes
from pathlib import Path

import numpy as np
from typing import Literal

from . import mat_index_comb

from .array_to_c_pointer_convert import arr_2d_to_c

from numpy.typing import NDArray
from numpy.ctypeslib import ndpointer

uint_c_type = np.uintp

ENUM_DCCA_of = Literal['all']

def p_dcca(input_data: NDArray[np.float64], 
           tws:  NDArray[np.int64] | NDArray[np.float64], 
           DCCA_of: np.ndarray | ENUM_DCCA_of ="all",
           axis = 0,
           P_DCCA_output_matrix: bool=False
        ) ->    tuple[
                NDArray[np.float64],    # DFA
                NDArray[np.float64],    # DCCA
                NDArray[np.float64]     # P_DCCA
                ]:

    """A function that calculates the 
       <span>&Rho;<sub>DCCA</sub></span>
        for a group of time series

        Args:
            input_data (NDArray[np.float64]): 2D array of times series integrated data.
            tws (NDArray[np.float64] | NDArray[np.float64]): 1D array of time scales.
            DCCA_of (np.ndarray | list | None, optional): _description_. Defaults to None.
            axis (int): axis of the input_data matrix that contains the values for each time series. Defalts to 0.
            P_DCCA_output_matrix (bool, optional): _description_. Defaults to False.

        Returns:
            <span>A tuple of 3 matrices:</span>
            DFA(NDArray[np.float64]):_description_. 
            DCCA(NDArray[np.float64]):_description_.
            <span>&Rho;<sub>DCCA</sub></span>(NDArray[np.float64]):_description_.
    """
    
    assert (tws[:-1] < tws[1:]).all() == True , ("""time window scales (tws) values must be in crescent order.""")


    # setting lib path

    sys_dict = {'Windows': 'zebendezig.dll', 'Linux': 'libzebendezig.so', 'Darwin': 'libzebendezig.dylib'}

    sys_info = platform.uname()

    lib_folder = Path(__file__).parent / 'zig_libs'

    lib_file = lib_folder / sys_dict[sys_info.system]

    zz = ctypes.cdll.LoadLibrary(lib_file)

    # outputs types
    c_2d_any_1d_uint = ndpointer(dtype = uint_c_type, ndim=1, flags='C')
    c_1d_double = ndpointer(dtype=np.double, ndim=1, flags='C')
    # getting the module file
    _p_dcca = zz.p_dcca
    # setting outputs
    _p_dcca.argtypes = [ 
                    c_1d_double, ctypes.c_size_t, ctypes.c_size_t, # input data
                    c_2d_any_1d_uint, ctypes.c_size_t, # tws
                    c_1d_double, # Time steps
                    c_2d_any_1d_uint, ctypes.c_size_t, # DCCA of
                    # outputs
                    c_2d_any_1d_uint, # DFA_arr
                    c_2d_any_1d_uint, # DCCA_arr
                    c_1d_double, ctypes.c_bool, # DCCA_arr
                    ]

    # Preparing input data    

    if not tws.flags.c_contiguous:
        tws = np.ascontiguousarray(tws)
    
    # preparing DCCA_ot array

    if type(DCCA_of) == str:
        if DCCA_of == "all":
            DCCA_of =  np.ascontiguousarray(mat_index_comb(input_data, axis=axis))
    # ensuring data compatibility
    DCCA_of = DCCA_of.astype(uint_c_type)
    c_DCCA_of = arr_2d_to_c(DCCA_of)

    
    # preparing tws array
    tws = tws.astype(uint_c_type)

    # preparing output array
    F_DFA_arr = np.ascontiguousarray(np.zeros(shape=(tws.shape[0], input_data.shape[axis]), dtype=input_data.dtype))
    c_DFA_arr = arr_2d_to_c(F_DFA_arr)
    DCCA_arr = np.ascontiguousarray(np.zeros(shape=(tws.shape[0], DCCA_of.shape[0]), dtype=input_data.dtype))
    c_DCCA_arr = arr_2d_to_c(DCCA_arr)

    # preparing data
    data_shape = input_data.shape

    if axis == 0:
        x_cnt = data_shape[0]
        x_len = data_shape[1]

        input_data = np.ascontiguousarray(input_data.flatten())

    if axis == 1:
        x_len = data_shape[0]
        x_cnt = data_shape[1]

        input_data = np.ascontiguousarray(np.asfortranarray(input_data).flatten(order="K"))

    time_steps = np.ascontiguousarray(np.arange(x_len, dtype=input_data.dtype))

    # preparing P_DCCA array
    # output array
    if P_DCCA_output_matrix == False:
        P_DCCA_shape = np.array([tws.shape[0], DCCA_of.shape[0]])
        P_DCCA_arr = np.ascontiguousarray(np.full(shape=(P_DCCA_shape.prod()) ,fill_value=np.nan, dtype=input_data.dtype))
    # output matrix
    elif P_DCCA_output_matrix == True:
        P_DCCA_arr = np.full(
            shape=(x_cnt, x_cnt, tws.shape[0]),fill_value=np.nan, dtype=input_data.dtype)
        # fill diagonal with ones
        P_DCCA_shape =  P_DCCA_arr.shape
        r = np.arange(x_cnt, dtype= np.int64)
        P_DCCA_arr[r,r, :] = 1
        del r
        # flatteing matrix
        P_DCCA_arr = np.ascontiguousarray(P_DCCA_arr.flatten())

    # calling functon  
    _p_dcca(input_data,  x_len, x_cnt, 
            tws, tws.size, 
            time_steps, 
            c_DCCA_of , DCCA_of.shape[0],
            # Outputs
            c_DFA_arr,
            c_DCCA_arr,
            P_DCCA_arr,
            # P_DCCA output matrix
            P_DCCA_output_matrix
            )

    # reshaping P_DCCA output
    P_DCCA_arr = P_DCCA_arr.reshape(P_DCCA_shape)

    return [F_DFA_arr, DCCA_arr, P_DCCA_arr]
