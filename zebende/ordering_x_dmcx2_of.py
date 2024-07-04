import numpy as np


def ordering_x_dmcx2_of(dmcx2_of):
    y_serie = dmcx2_of[:, 0]
    x_series = np.sort(dmcx2_of[:, 1:], axis=1)
    out = np.c_[y_serie, x_series]
    return out
