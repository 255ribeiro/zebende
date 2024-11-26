import numpy as np
from numpy.typing import NDArray


# time window array
def time_windows(n_min: int, n_max: int, exp_fac:float = 8.0, ensure_max:bool = True) -> NDArray[np.int64]:

    """Generate an ordered array of time sclaes for the Detrended calculations.
    Args:
        n_min (int): Minimun value for a time window (box).
        n_max (int): Maximun value for a time window (box).
        exp_fac (float, optional): Exponencial factor to calculate successive time windows in the array. Defaults to 8.0.
        ensure_max (bool, optional): Ensure that the Maximun value is. Defaults to True.
    Returns:
        NDArray[np.int64]: Time scales array.
    """
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
