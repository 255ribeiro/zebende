import numpy as np
import numpy.typing as npt


# time window array
def time_windows(n_min: np.int32, n_max: int, exp_fac: np.float64 = 8.0) -> npt.NDArray[np.float64]:
    n = n_min
    tmp = []
    ir = 0

    while n <= n_max:
        tmp.append(n)
        ir = ir + 1
        n = int((n_min + ir) * np.power(np.power(2, 1.0 / exp_fac), ir))
    if tmp[-1] < n_max:
        tmp.append(n_max)

    return np.array(tmp)
