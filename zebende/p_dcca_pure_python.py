import numpy as np
from typing import Literal

from numpy.typing import NDArray

from . import detrended_series, mat_index_comb

from .p_dcca_simple_output import p_dcca_simple_output
from .p_dcca_matrix_output import p_dcca_matrix_output

ENUM_DCCA_of = Literal['all']

# P_DCCA calculator
def p_dcca(
    data: NDArray[np.float64], tws: NDArray[np.int64] | NDArray[np.float64], DCCA_of: np.ndarray | ENUM_DCCA_of ="all", P_DCCA_output_format="simple"
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    # setting time_steps
    time_steps = np.arange(data.shape[0])
    if type(DCCA_of) == str:
        if DCCA_of == "all":
            DCCA_of = mat_index_comb(data, axis=1)
    # Global outputs for DFA and DCCA
    F_DFA_arr = np.zeros(shape=(tws.shape[0], data.shape[1]), dtype=data.dtype)
    DCCA_arr = np.zeros(shape=(tws.shape[0], DCCA_of.shape[0]), dtype=data.dtype)
    # P_DCCA global output format and function
    # simple output
    if P_DCCA_output_format == "simple":
        # output array
        P_DCCA_arr = np.full(shape=(tws.shape[0], DCCA_of.shape[0]) ,fill_value=np.nan, dtype=data.dtype)
        # output function
        P_DCCA_output_funtion = p_dcca_simple_output
    # matrix output
    elif P_DCCA_output_format == "matrix":
        # output matrix
        P_DCCA_arr = np.full(
            shape=(DCCA_of.max() + 1, DCCA_of.max() + 1, tws.shape[0]),fill_value=np.nan, dtype=data.dtype
        )
        # fill diagonal with ones
        r = np.arange(DCCA_of.max() + 1)
        P_DCCA_arr[r,r, :] = 1
        del r
        # output function
        P_DCCA_output_function = p_dcca_matrix_output
    # for time scales in n
    for n_index, n  in enumerate(tws):
        # in time scale (n) accumulators
        f2dfa_n = np.full(shape=(data.shape[0] - n, data.shape[1]),fill_value=np.nan, dtype=data.dtype)
        dcca_n = np.full(shape=(data.shape[0] - n, DCCA_of.shape[0]), fill_value=np.nan, dtype=data.dtype)
        detrended_mat = np.full(shape=(n + 1, data.shape[1]), fill_value=np.nan, dtype=data.dtype)

        for i in range(data.shape[0] - n): # for each box

            detrended_series( # inputs
                time_steps[i : i + (n + 1)],  # arr_x
                data[i : i + (n + 1), :],  # mat_y
                detrended_mat,  # output
            )

            f2dfa_n[i] = np.power(detrended_mat, 2).mean(axis=0)

            for j, pair in enumerate(DCCA_of):

                dcca_n[i, j] = (
                    detrended_mat[:, pair[0]] * detrended_mat[:, pair[1]]
                ).mean(axis=0)

        F_DFA_arr[n_index, :] = np.sqrt(f2dfa_n.mean(axis=0))

        DCCA_arr[n_index, :] = dcca_n.mean(axis=0)

        # calculation of P_DCCA
        P_DCCA_output_function(n_index, DCCA_of, F_DFA_arr, DCCA_arr, P_DCCA_arr)

    return F_DFA_arr, DCCA_arr, P_DCCA_arr
