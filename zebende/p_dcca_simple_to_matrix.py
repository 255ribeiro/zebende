import numpy as np
from . import mat_index_comb
from numpy.typing import NDArray

def p_dcca_simple_to_matrix(simple_pdcca:NDArray[np.float64], DCCA_of:NDArray[np.int64] | None = None ) ->NDArray[np.float64]:

    """Convert a simple (table) of  <span>&Rho;<sub>DCCA</sub></span> values to the matrix output.

    Args:
        simple_pdcca:NDArray[np.float64]: Table of <span>&Rho;<sub>DCCA</sub></span> outputs.
        DCCA_of (np.ndarray | None, optional): The DCCA_of used for the <span>&Rho;<sub>DCCA</sub></span> calculations.
        If None, the function will assume it was created with the Default option in the p_dcca function. Defaults to None.

    Returns:
        NDArray[np.float64]: 3D matrix with the values of the <span>&Rho;<sub>DCCA</sub></span> of each time scale as a level.
    """
    
    if DCCA_of is None:
        series_count = int(np.round((1 + np.sqrt(1 + 8 *simple_pdcca.shape[1]))/2, 0))
    
        DCCA_of = mat_index_comb(series_count)
    else:
        series_count = DCCA_of.flatten().max()+1

    tws_count = simple_pdcca.shape[0]

    shape = (series_count, series_count, tws_count)
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
