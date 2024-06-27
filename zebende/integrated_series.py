import numpy as np
import numpy.typing as npt


# integrates series
def integrated_series(mat_series):
    out = (mat_series - mat_series.mean(axis=0)).cumsum(axis=0)
    return out
