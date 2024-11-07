from itertools import combinations

import numpy as np

from numpy.typing import NDArray

def mat_index_comb(mat:NDArray[np.float64] | int, axis:int=0) -> NDArray[np.int64] | None:
    if type(mat) == np.ndarray:
        return np.array(list(combinations(range(mat.shape[axis]), 2)))
    elif type(mat) == int:
        return np.array(list(combinations(range(mat), 2)))
    else:
        print('mat of type {} not supported'.format(type(mat)))
