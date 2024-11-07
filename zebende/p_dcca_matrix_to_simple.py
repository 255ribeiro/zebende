import numpy as np

from . import mat_index_comb

from numpy.typing import NDArray

def p_dcca_matrix_to_simple(matrix_pdcca:NDArray[np.float64]) ->NDArray[np.float64]:
    series_count = matrix_pdcca.shape[0]
    tws_count = matrix_pdcca.shape[2]
    DCCA_of = mat_index_comb(series_count)
    shape = (tws_count, DCCA_of.shape[0])
   

    P_DCCA_arr = np.full(shape=shape,fill_value=np.nan, dtype=matrix_pdcca.dtype)

    for n_index in range(tws_count):
        for i, pair in enumerate(DCCA_of):
            P_DCCA_arr[n_index, i] = matrix_pdcca[pair[0], pair[1], n_index] 

    return P_DCCA_arr
