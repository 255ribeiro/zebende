
import numpy as np
from numba import prange
from numpy.typing import NDArray

from . import detrended_series, mat_index_comb

from .p_dcca_simple_output import p_dcca_simple_output
from .p_dcca_matrix_output import p_dcca_matrix_output


# P_DCCA calculator
def p_dcca(data, tws, time_steps=None, DCCA_of='all', P_DCCA_output_format = 'simple')-> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:

    # setting time_steps in None is passed
    if time_steps == None:
        time_steps = np.arange(data.shape[0])

    if type(DCCA_of) == str:
        if DCCA_of == 'all':
            DCCA_of = mat_index_comb(data, axis=1)
    # Global outputs
    F_DFA_arr = np.empty(shape=(tws.shape[0], data.shape[1]), dtype=data.dtype)
    DCCA_arr = np.empty(shape=(tws.shape[0], DCCA_of.shape[0]), dtype=data.dtype)

    # P_DCCA output formats and functions

    # simple output
    if P_DCCA_output_format == 'simple':
        #output array
        P_DCCA_arr = np.empty(shape=(tws.shape[0], DCCA_of.shape[0]), dtype=data.dtype)
        # output function
        P_DCCA_output_funtion = p_dcca_simple_output
    # matrix output
    elif P_DCCA_output_format == 'matrix':
        # output matrix
        P_DCCA_arr = np.ones(shape=(DCCA_of.max()+1, DCCA_of.max()+1, tws.shape[0]), dtype=data.dtype)
        # output function
        P_DCCA_output_funtion = p_dcca_matrix_output


    # for time scales in n
    for n_index in range(len(tws)):

        n = tws[n_index]

        # in time scale (n) accumulators
        dcca_n = np.zeros(
            shape=(data.shape[0] - n, DCCA_of.shape[0]), dtype=data.dtype
        )

        f2dfa_n = np.empty(
            shape=(data.shape[0] - n, data.shape[1]), dtype=data.dtype
        )

        # Operações dentro das caixas (sobrepostas)

        for i in range(data.shape[0] - n):

            # 2-- dividir o sinal em N-n janelas temporais de comprimento n
            # janela temporal
            time_window = data[i : i + (n + 1), :]

            time_steps_window = time_steps[i : i + (n + 1)]

            # 3-- Ajustar uma curva de tedência

            # geralente polinômio de primerio grau

            detrended_mat = detrended_series(time_steps_window, time_window)

            f2dfa_n[i] = np.power(detrended_mat, 2).mean(axis=0)

            for j, comb in enumerate(DCCA_of):

                dcca_n[i, j] = (detrended_mat[:, comb[0]] * detrended_mat[:, comb[1]]).mean(axis=0)

        # 5--para cada escala temporal
        # salvar valor de cada escala temporal

        F_DFA_arr[n_index, :] = np.sqrt(f2dfa_n.mean(axis=0))

        DCCA_arr[n_index, :] = dcca_n.mean(axis=0)

        # calculation of P_DCCA
        P_DCCA_output_funtion(n_index, DCCA_of, F_DFA_arr, DCCA_arr, P_DCCA_arr)

    return F_DFA_arr, DCCA_arr, P_DCCA_arr
