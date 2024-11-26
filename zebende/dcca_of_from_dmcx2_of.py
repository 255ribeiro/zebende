from itertools import combinations
from numpy.typing import NDArray

import numpy as np


def dcca_of_from_dmcx2_of(dmc_of:NDArray[np.int64]) -> NDArray[np.int64]:
    """_summary_

    Args:
        dmc_of (NDArray[np.int64]): _description_

    Returns:
        NDArray[np.int64]: _description_
    """    
    if type(dmc_of) == np.ndarray:
        dmc_of = dmc_of.tolist()
    out = []
    for i in range(len(dmc_of)):
        temp = list(combinations(dmc_of[i], 2))
        for j in temp:
            j = list(j)

            j.sort()

            if j not in out:
                out.append(j)
    return np.array(out)
