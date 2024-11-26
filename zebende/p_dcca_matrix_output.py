import numpy as np

def p_dcca_matrix_output(n_index:int, DCCA_of, F_DFA_arr, DCCA_arr, P_DCCA_arr):  

    for i in range(DCCA_of.shape[0]):
            p_dcca = DCCA_arr[n_index, i] / (F_DFA_arr[n_index, DCCA_of[i][0]] * F_DFA_arr[n_index, DCCA_of[i][1]])
            P_DCCA_arr[DCCA_of[i, 0], DCCA_of[i, 1], n_index] = p_dcca
            P_DCCA_arr[DCCA_of[i, 1], DCCA_of[i, 0], n_index] = p_dcca
