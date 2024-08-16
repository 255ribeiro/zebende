import numpy as np

from numpy.typing import NDArray

from . import detrended_series, mat_index_comb, detrended_fit_series

from .p_dcca_simple_output import p_dcca_simple_output
from .p_dcca_matrix_output import p_dcca_matrix_output


# P_DCCA calculator
def dfa(data, tws, time_steps=None ) -> NDArray[np.float64]:
        # setting time_steps in None is passed
    if time_steps == None:
        time_steps = np.arange(data.shape[0])


    # Global outputs
    F_DFA_arr = np.zeros(shape=(tws.shape[0], data.shape[1]), dtype=data.dtype)


    # for time scales in n
    for n_index in range(len(tws)):

        n = tws[n_index]

        # in time scale (n) accumulators


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



        # 5--para cada escala temporal
        # salvar valor de cada escala temporal

        F_DFA_arr[n_index, :] = np.sqrt(f2dfa_n.mean(axis=0))



    return F_DFA_arr
