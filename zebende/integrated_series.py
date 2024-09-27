import numpy as np

# integrates series
def integrated_series(mat_series):
    out = (mat_series - mat_series.mean(axis=0)).cumsum(axis=0)
    return out
