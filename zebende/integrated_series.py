import numpy as np
from numpy.typing import NDArray

# integrates series
def integrated_series(mat_series:NDArray[np.float64], axis:int = 0, return_mean:bool = False) -> NDArray[np.float64]:
    """Returns a matrix of integrates series from a matrix of time series.
        The integrates series is a cumulative sum of the values of the series subtracted by the mean.

    Args:
        mat_series (NDArray[np.float64]): Matrix of time series with one serie per column.
        axis (int): axis of the input_data matrix that contains the values for each time series. Defalts to 0.
        return_mean (bool): If an array containing the means of each series should be returned. Defaults to False.

    Returns:
        NDArray[np.float64]: Matrix of integrated time series with one integrated time series per column.
        NDArray[np.float64]: Returned only if return_mean == True. Array with the mean value of each feature.
    """

    if axis == 0:
        series_means = mat_series.mean(axis=1)
        out = (mat_series - series_means.reshape(series_means.size, 1)).cumsum(axis=1)

    elif axis == 1:
        series_means = mat_series.mean(axis=0)
        out = (mat_series - series_means).cumsum(axis=0)
    
    if return_mean == True:
        return out,  series_means
    else:
        return out
