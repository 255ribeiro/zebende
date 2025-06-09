import numpy as np
from numpy.typing import NDArray


def ordering_x_dmcx2_of(dmcx2_of:NDArray[np.int64])-> NDArray[np.int64]:
    """Ordering the independent indexes of the dmcx2_of NDArray

    Args:
        dmcx2_of (NDArray[np.int64]): dmcx2_of NDArray

    Returns:
        NDArray[np.int64]: dmcx2_of NDArray with the independent variables ordered.
    """
    
    if   type(dmcx2_of) == np.ndarray:
        if dmcx2_of.ndim == 2:
            y_serie = dmcx2_of[:, 0]
            x_series = np.sort(dmcx2_of[:, 1:], axis=1)
            out = np.c_[y_serie, x_series]

        elif dmcx2_of.ndim == 1:
            y_serie = dmcx2_of[0:1]
            x_series = np.sort(dmcx2_of[1:])
            out =  np.concatenate((y_serie, x_series))
        else:
            print("dmcx2_of must be a 2D or 1D Numpy array")

    return out
