import numpy as np

def p_dcca_simple_output(n_index, DCCA_of, F_DFA_arr, DCCA_arr, P_DCCA_arr):

    for i, dcca_pair in enumerate(DCCA_of):

                P_DCCA_arr[n_index, i] = DCCA_arr[n_index, i] / (F_DFA_arr[n_index, dcca_pair[0]] * F_DFA_arr[n_index, dcca_pair[1]])
    return P_DCCA_arr
    