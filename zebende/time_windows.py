import numpy as np
from numpy.typing import NDArray


# time window array
def time_windows(n_min: np.int32, n_max: int, exp_fac: np.float64 = 8.0, ensure_max = True) -> NDArray[np.int64]:
    n = n_min
    tmp = []
    ir = 0
    while n <= n_max:
        tmp.append(n)
        ir = ir + 1
        n = int((n_min + ir) * np.power(np.power(2, 1.0 / exp_fac), ir))
    if ((ensure_max == True) and (tmp[-1] < n_max)):
        tmp.append(n_max)

    return np.array(tmp, dtype=np.int64)
