from typing import (Literal, Any, Union)

import numpy as np
from numpy.typing import (NDArray, DTypeLike)


from . import (
    dcca_of_from_dmcx2_of,
    dmc_of_all_as_y,
    dmcx2_from_p_dcca_matrix,
    p_dcca_pure_python,
)

ENUM_DMCx2_of = Literal['all-full', 'first-full'] 


def dmcx2_pure_python(input_data: NDArray[np.float64], 
          tws: NDArray[np.int64] | NDArray[np.float64], 
          dmcx2_of: NDArray[np.float64] | list | ENUM_DMCx2_of = 'all-full'
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
            dmcx2_of (NDArray[np.float64] | list | ENUM_DMCx2_of, optional): _description_. Defaults to 'all-full'.
            DCCA_of (np.ndarray | list | None, optional): _description_. Defaults to None.

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
            dmcx2_of = np.array( [np.arange(input_data.shape[1])])
        # creating ndarray of y and x values for DMCx2 calculations
        elif dmcx2_of == 'all-full':
            dmcx2_of = dmc_of_all_as_y(input_data)
    
    test_dmcx2_of = dmcx2_of[:,1:]
    assert (test_dmcx2_of[:,:-1] < test_dmcx2_of[:,1:]).all() == True , ("""
Dmcx2 x values out of order: use zebende.ordering_x_dmcx2_of(dmcx2_of) to fix it before passing the dmcx2_of value to zebende.dmcx2() function""")
    del test_dmcx2_of

    # creating ndarray for P_DCCA calculations based on the DMCx2 array
    DCCA_of = dcca_of_from_dmcx2_of(dmcx2_of)


    # P_DCCA calculations
    F_DFA_arr, DCCA_arr, P_DCCA_arr = p_dcca_pure_python(input_data=input_data, tws=tws, DCCA_of=DCCA_of,  P_DCCA_output_matrix = True)

    # DMCx2 output matrix
    DMCx2_arr = np.full(shape=(tws.shape[0], dmcx2_of.shape[0]), fill_value=np.nan, dtype=input_data.dtype)

    DMCx2_arr = dmcx2_from_p_dcca_matrix(P_DCCA_arr, tws, dmcx2_of)

    return F_DFA_arr, DCCA_arr, P_DCCA_arr, DMCx2_arr
