from typing import (Literal, Any, Union)

import numpy as np
from numpy.typing import (NDArray, DTypeLike)


from . import (
    dcca_of_from_dmcx2_of,
    dmc_of_all_as_y,
    dmcx2_from_p_dcca_matrix,
    p_dcca,
)

ENUM_DMCx2_of = Literal['all-full', 'first-full'] 


def dmcx2(input_data: NDArray[np.float64], 
          tws: NDArray[np.int64] | NDArray[np.float64], 
          dmcx2_of: NDArray[np.float64] | Literal['all-full', 'first-full'] = 'all-full',
           axis: int = 0
           )-> tuple[
                                                            NDArray[np.float64],
                                                            NDArray[np.float64],
                                                            NDArray[np.float64],
                                                            NDArray[np.float64]
                                                             ]:

    """
        A function that calculates the <span>DMC<sub>x</sub><sup>2</sup></span> for a group of time series

        Args:
            input_data (NDArray[np.float64]): _description_.
            tws (NDArray[np.int64] | NDArray[np.float64]): _description_.
            dmcx2_of (NDArray[np.float64] | Literal['all-full', 'first-full'], optional): _description_. Defaults to 'all-full'.
            DCCA_of (np.ndarray | list | None, optional): _description_. Defaults to None.
            axis (int): axis of the input_data matrix that contains the values for each time series. Defalts to 0.

        Returns:
            <span>A tuple of 4 matrices:</span><br>
            DFA(NDArray[np.float64]):_description_,<br>
            DCCA(NDArray[np.float64]):_description_,<br>
            <span>&Rho;<sub>DCCA</sub></span>(NDArray[np.float64]):_description_,<br>
            <span>DMC<sub>x</sub><sup>2</sup></span>(NDArray[np.float64]):_description_.<br>

    """
 
    
    if type(dmcx2_of) == str:

        # creating ndarray of y and x values for DMCx2 calculations
        if dmcx2_of == 'first-full':
            dmcx2_of = np.array( [np.arange(input_data.shape[axis])])
        # creating ndarray of y and x values for DMCx2 calculations
        elif dmcx2_of == 'all-full':
            dmcx2_of = dmc_of_all_as_y(input_data, axis=axis)
    

    # creating ndarray for P_DCCA calculations based on the DMCx2 array
    DCCA_of = dcca_of_from_dmcx2_of(dmcx2_of)

    # P_DCCA calculations
    F_DFA_arr, DCCA_arr, P_DCCA_arr = p_dcca(input_data=input_data, tws=tws, DCCA_of=DCCA_of,  P_DCCA_output_matrix = True)

    DMCx2_arr = dmcx2_from_p_dcca_matrix(P_DCCA_arr, dmcx2_of)

    return F_DFA_arr, DCCA_arr, P_DCCA_arr, DMCx2_arr
