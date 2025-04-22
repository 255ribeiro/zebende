import numpy as np
from . import ordering_x_dmcx2_of
from numpy.typing import NDArray


def dmc_of_all_as_y(data:NDArray[np.float64],
                    axis:int = 0,
                    ) -> NDArray[np.float64]:
    """Generates a 2D array for  <span>DMC<sub>x</sub><sup>2</sup></span> calculations

    Args:
        data (NDArray[np.float64]): 2D array of data.

    Returns:
        NDArray[np.float64]: 2D array where in each row, the first element is the dependent variable,
        and the other valus are the intependent ones.
        In each line one of the series is choosen as the dependent variable,
        so the number of lines correspond to the number of columns in the data array.
    """    
    aux = list(range(data.shape[axis]))
    dmc_list = [aux]
    for i in range(data.shape[axis]-1):
        aux = aux[1:] + aux[:1]
        dmc_list.append(aux)

    dmc_list = np.array(dmc_list)

    dmc_list = ordering_x_dmcx2_of(dmc_list)

    return dmc_list

