import numpy as np

from numpy.typing import NDArray

from . import detrended_series

# P_DCCA calculator
def dfa_pure_python(data: NDArray[np.float64], tws: NDArray[np.int64]) -> NDArray[np.float64]:
    
    """Calculates the Detrended Fluctuation analysis (DFA) for one or more time series. 

    Args:
        data (NDArray[np.float64]): data matrix with one series per column.
        tws (NDArray[np.int64]): 1D array of time windows.

    Returns:
        NDArray[np.float64]: DFA table.
    """ 
        # setting time_steps in None is passed

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
