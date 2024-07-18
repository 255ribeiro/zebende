import numpy as np
from numba import prange
from numpy.typing import NDArray

from . import detrended_series, mat_index_comb, detrended_fit_series

from .p_dcca_simple_output import p_dcca_simple_output
from .p_dcca_matrix_output import p_dcca_matrix_output


# P_DCCA calculator
def p_dcca(
    data, tws, time_steps=None, DCCA_of="all", P_DCCA_output_format="simple"
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:

    # setting time_steps in None is passed
    if time_steps == None:
        time_steps = np.arange(data.shape[0])

    if type(DCCA_of) == str:
        if DCCA_of == "all":
            DCCA_of = mat_index_comb(data, axis=1)
    # Global outputs
    F_DFA_arr = np.empty(shape=(tws.shape[0], data.shape[1]), dtype=data.dtype)
    DCCA_arr = np.empty(shape=(tws.shape[0], DCCA_of.shape[0]), dtype=data.dtype)

    # P_DCCA output formats and functions

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
        np.fill_diagonal(P_DCCA_arr, val=1)

        # output function
        P_DCCA_output_funtion = p_dcca_matrix_output

    # for time scales in n
    for n_index in range(len(tws)):

        n = tws[n_index]

        # in time scale (n) accumulators
        dcca_n = np.full(shape=(data.shape[0] - n, DCCA_of.shape[0]), fill_value=np.nan, dtype=data.dtype)

        f2dfa_n = np.full(shape=(data.shape[0] - n, data.shape[1]),fill_value=np.nan, dtype=data.dtype)

        detrended_mat = np.full(shape=(n + 1, data.shape[1]), fill_value=np.nan, dtype=data.dtype)
    
        # Operações dentro das caixas (sobrepostas)

        for i in range(data.shape[0] - n):

            # 2-- dividir o sinal em N-n janelas temporais de comprimento n
            # janela temporal

            # 3-- Ajustar uma curva de tedência

            # geralente polinômio de primerio grau

            detrended_series( # inputs
                time_steps[i : i + (n + 1)],  # arr_x
                data[i : i + (n + 1), :],  # mat_y
                detrended_mat,  # output
            )

            f2dfa_n[i] = np.power(detrended_mat, 2).mean(axis=0)

            for j, comb in enumerate(DCCA_of):

                dcca_n[i, j] = (
                    detrended_mat[:, comb[0]] * detrended_mat[:, comb[1]]
                ).mean(axis=0)

        # 5--para cada escala temporal
        # salvar valor de cada escala temporal

        F_DFA_arr[n_index, :] = np.sqrt(f2dfa_n.mean(axis=0))

        DCCA_arr[n_index, :] = dcca_n.mean(axis=0)

        # calculation of P_DCCA
        P_DCCA_output_funtion(n_index, DCCA_of, F_DFA_arr, DCCA_arr, P_DCCA_arr)

    return F_DFA_arr, DCCA_arr, P_DCCA_arr
