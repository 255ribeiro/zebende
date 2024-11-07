import numpy as np
from . import mat_index_comb
from numpy.typing import NDArray

def p_dcca_simple_to_matrix(simple_pdcca:NDArray[np.float64]) ->NDArray[np.float64]:

    series_count = int(np.round((1 + np.sqrt(1 + 8 *simple_pdcca.shape[1]))/2, 0))
    tws_count = simple_pdcca.shape[0]
    DCCA_of = mat_index_comb(series_count)
    shape = (series_count, series_count, tws_count)
    print(shape)
    P_DCCA_arr = np.full(shape=shape,fill_value=np.nan, dtype=simple_pdcca.dtype)
    # fill diagonal with ones
    r = np.arange(series_count, dtype= np.int64)
    P_DCCA_arr[r,r, :] = 1
    del r

    for n_index in range(tws_count):
        for i, pair in enumerate(DCCA_of):
            P_DCCA_arr[pair[0], pair[1], n_index] = simple_pdcca[n_index, i]
            P_DCCA_arr[pair[1], pair[0], n_index] = simple_pdcca[n_index, i]

    return P_DCCA_arr
