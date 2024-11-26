import numpy as np
from numpy.typing import NDArray

# integrates series
def integrated_series(mat_series:NDArray[np.float64]) -> NDArray[np.float64]:
    """Returns a matrix of integrates series from a matrix of time series.

    Args:
        mat_series (NDArray[np.float64]): Matrix of time series with one serie per column.

    Returns:
        NDArray[np.float64]: Matrix of integrated time series with one integrated time series per column.
    """    
    out = (mat_series - mat_series.mean(axis=0)).cumsum(axis=0)
    return out
