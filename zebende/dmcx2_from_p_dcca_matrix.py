from typing import (Literal, Any, Union)

import numpy as np
from numpy.typing import (NDArray, DTypeLike)


ENUM_DMCx2_of = Literal['all-full', 'first-full'] 


# DMCx2 output matrix
def dmcx2_from_p_dcca_matrix(P_DCCA_arr: NDArray[np.float64], tws: NDArray[np.float64],  dmcx2_of: NDArray[np.float64]) -> NDArray[np.float64]:

    DMCx2_arr = np.full(shape=(tws.shape[0], dmcx2_of.shape[0]), fill_value=np.nan, dtype=P_DCCA_arr.dtype)

    for n_index in range(len(tws)):

        for j in range(dmcx2_of.shape[0]):

            y_indexes = dmcx2_of[j, 0:1]
            x_indexes = dmcx2_of[j, 1:]


            mat_x = P_DCCA_arr[np.ix_(x_indexes, x_indexes)][:, :, n_index]
            vec_y = P_DCCA_arr[np.ix_(x_indexes, y_indexes)][:, :, n_index]
            # print(mat_x, vec_y)
            # dmcx2 calculation
            DMCx2_arr[n_index, j] = vec_y.T @ np.linalg.inv(mat_x) @ vec_y
    return DMCx2_arr
