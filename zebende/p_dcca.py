import numpy as np
import numpy.typing as npt


# P_DCCA calculator
def p_dcca(int_series_mat , tws, time_stamps=None, DCCA_of=None):
    # Global outputs
    F_DFA_arr = np.empty(shape=(tws.shape[0], int_series_mat.shape[1]), dtype=float)
    DCCA_arr = np.empty(shape=(tws.shape[0], DCCA_of.shape[0]), dtype=float)
    P_DCCA_arr = np.empty(shape=(tws.shape[0], DCCA_of.shape[0]), dtype=float)

    # for time scales in n
    for n_index in range(len(tws)):

        n = tws[n_index]

        # in time scale (n) accumulators
        dcca_n = np.zeros(
            shape=(int_series_mat.shape[0] - n, DCCA_of.shape[0]), dtype=np.float64
        )

        f2dfa_n = np.empty(
            shape=(int_series_mat.shape[0] - n, int_series_mat.shape[1]), dtype=np.float64
        )

        # Operações dentro das caixas (sobrepostas)

        for i in range(int_series_mat.shape[0] - n):

            # 2-- dividir o sinal em N-n janelas temporais de comprimento n
            # janela temporal
            time_window = int_series_mat[i : i + (n + 1), :]

            time_stamps_window = time_stamps[i : i + (n + 1)]

            # 3-- Ajustar uma curva de tedência

            # geralente polinômio de primerio grau

            detrended_mat = detrended_serie(time_stamps_window, time_window)

            f2dfa_n[i] = np.power(detrended_mat, 2).mean(axis=0)

            for j, comb in enumerate(DCCA_of):

                dcca_n[i, j] = (detrended_mat[:, comb[0]] * detrended_mat[:, comb[1]]).mean(axis=0)

        # 5--para cada escala temporal
        # salvar valor de cada escala temporal

        F_DFA_arr[n_index, :] = np.sqrt(f2dfa_n.mean(axis=0))

        DCCA_arr[n_index, :] = dcca_n.mean(axis=0)

        for j, l in enumerate(DCCA_of):

            P_DCCA_arr[n_index, j] = np.round(DCCA_arr[n_index, j] / (F_DFA_arr[n_index, l[0]] * F_DFA_arr[n_index, l[1]]), 5)

    return F_DFA_arr, DCCA_arr, P_DCCA_arr
