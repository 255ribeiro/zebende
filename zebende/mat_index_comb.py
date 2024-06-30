from itertools import combinations

import numpy as np


def mat_index_comb(mat, axis=0):

    return np.array(list(combinations(range(mat.shape[axis]), 2)))
