from itertools import combinations

import numpy as np

from numpy.typing import NDArray

def mat_index_comb(mat:NDArray[np.float64] | int, axis:int = 0) -> NDArray[np.uintp] | None:

    """_summary_

    Args:
        mat (NDArray[np.float64], int ): Matrix of time series with one serie per column.
        axis (int , optional):

    Returns:
        (NDArray[np.uintp], None): _description_
    """
    if type(mat) == np.ndarray:
        return np.array(list(combinations(range(mat.shape[axis]), 2)), dtype=np.uintp)
    elif type(mat) == int:
        return np.array(list(combinations(range(mat), 2)), dtype=np.uintp)
    else:
        print('mat of type {} not supported'.format(type(mat)))
        return None
