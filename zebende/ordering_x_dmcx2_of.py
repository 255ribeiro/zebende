import numpy as np
from numpy.typing import NDArray


def ordering_x_dmcx2_of(dmcx2_of:NDArray[np.int64])-> NDArray[np.int64]:
    """_summary_

    Args:
        dmcx2_of (NDArray[np.int64]): _description_

    Returns:
        NDArray[np.int64]: _description_
    """    
    y_serie = dmcx2_of[:, 0]
    x_series = np.sort(dmcx2_of[:, 1:], axis=1)
    out = np.c_[y_serie, x_series]
    return out
