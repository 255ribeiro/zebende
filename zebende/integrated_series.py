import numpy as np
from numpy.typing import NDArray

# integrates series
def integrated_series(mat_series:NDArray[np.float64], axis:int = 0) -> NDArray[np.float64]:
    """Returns a matrix of integrates series from a matrix of time series.
        The integrates series is a cumulative sum of the values of the series subtracted by the mean.

    Args:
        mat_series (NDArray[np.float64]): Matrix of time series with one serie per column.
        axis (int): axis of the input_data matrix that contains the values for each time series. Defalts to 0.

    Returns:
        NDArray[np.float64]: Matrix of integrated time series with one integrated time series per column.
    """

    if axis == 0:
        out = (mat_series - mat_series.mean(axis=1)).cumsum(axis=1)

    elif axis == 1: 
        out = (mat_series - mat_series.mean(axis=0)).cumsum(axis=0)
    
    return out
