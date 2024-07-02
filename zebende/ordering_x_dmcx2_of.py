import numpy as np


def ordering_x_dmcx2_of(dmcx2_of):
    y_series = dmcx2_of[0, :]
    x_serie = dmcx2_of[1:, :]
    out = np.concat(np.sort(x_serie, axis=1), y_series, axis=1)
    return out