from typing import (Literal)

import numpy as np
from numpy.typing import (NDArray)


def dmcx2_from_p_dcca_matrix(P_DCCA_arr: NDArray[np.float64],  dmcx2_of: NDArray[np.float64]) -> NDArray[np.float64]:
    """_summary_

    Args:
        P_DCCA_arr (NDArray[np.float64]): _description_
        dmcx2_of (NDArray[np.float64]): _description_

    Returns:
        NDArray[np.float64]: _description_
    """    
    # DMCx2 output matrix
    DMCx2_arr = np.full(shape=(P_DCCA_arr.shape[2], dmcx2_of.shape[0]), fill_value=np.nan, dtype=P_DCCA_arr.dtype)

    for n_index in range(P_DCCA_arr.shape[2]):
        P_DCCA_arr_2D = P_DCCA_arr[:, :, n_index]
        for j,  dmcx2_of_1D in enumerate(dmcx2_of):
            DMCx2_arr[n_index, j] = dmcx2_from_p_dcca_matrix_2d(P_DCCA_arr_2D,  dmcx2_of_1D)
    return DMCx2_arr


def dmcx2_from_p_dcca_matrix_2d(P_DCCA_arr_2D: NDArray[np.float64],  dmcx2_of_1D: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculates the 

    Args:
        P_DCCA_arr_2D (NDArray[np.float64]): 2D array of <span>&Rho;<sub>DCCA</sub></span> results 
        dmcx2_of__1D (NDArray[np.float64]): 1D array of indexes

    Returns:
        <span>DMC<sub>x</sub><sup>2</sup></span>([float64]): _description_
    """    
    y_indexes = dmcx2_of_1D[0:1]
    x_indexes = dmcx2_of_1D[1:]

    mat_x = P_DCCA_arr_2D[np.ix_(x_indexes, x_indexes)]
    vec_y = P_DCCA_arr_2D[np.ix_(x_indexes, y_indexes)]

    return vec_y.T @ np.linalg.inv(mat_x) @ vec_y