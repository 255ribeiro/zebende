from typing import Literal

import numpy as np
from numpy.typing import NDArray


from . import (
    dcca_of_from_dmcx2_of,
    dmc_of_all_as_y,
    p_dcca,
)

ENUM_DMCx2_of = Literal['all-full', 'first-full']


def dmcx2(input_data: NDArray[np.float64], tws: NDArray[np.float64], dmcx2_of: NDArray[np.float64] | list | ENUM_DMCx2_of = 'all-full', time_steps: np.ndarray | None = None, DCCA_of: np.ndarray | list | None = None):
    
    if type(dmcx2_of) == str:

        # creating ndarray of y and x values for DMCx2 calculations
        if dmcx2_of == 'first-full':
            dmcx2_of =np.array( [np.arange(input_data.shape[1])])
        # creating ndarray of y and x values for DMCx2 calculations
        elif dmcx2_of == 'all-full':
            dmcx2_of = dmc_of_all_as_y(input_data)
    
    test_dmcx2_of = dmcx2_of[:,1:]
    assert (test_dmcx2_of[:,:-1] < test_dmcx2_of[:,1:]).all() == True , ("""
Dmcx2 x values out of order: use zebende.ordering_x_dmcx2_of(dmcx2_of) to fix it before passing the dmcx2_of value to zebende.dmcx2() function""")
    del test_dmcx2_of

    # creating ndarray for P_DCCA calculations based on the DMCx2 array
    if DCCA_of == None:
        DCCA_of = dcca_of_from_dmcx2_of(dmcx2_of)
    # P_DCCA calculations
    F_DFA_arr, DCCA_arr, P_DCCA_arr = p_dcca(input_data=input_data, tws=tws, time_steps=time_steps, DCCA_of=DCCA_of,  P_DCCA_output_matrix = True)

    # DMCx2 output matrix
    DMCx2_arr = np.empty(shape=(tws.shape[0], dmcx2_of.shape[0]), dtype=input_data.dtype)

    for n_index in range(len(tws)):

        for j in range(dmcx2_of.shape[0]):

            y_indexes = dmcx2_of[j, 0:1]
            x_indexes = dmcx2_of[j, 1:]

            mat_x = P_DCCA_arr[np.ix_(x_indexes, x_indexes)][:, :, n_index]
            vec_y = P_DCCA_arr[np.ix_(x_indexes, y_indexes)][:, :, n_index]
            # print(mat_x, vec_y)
            # dmcx2 calculation
            DMCx2_arr[n_index, j] = vec_y.T @ np.linalg.inv(mat_x) @ vec_y

    return F_DFA_arr, DCCA_arr, P_DCCA_arr, DMCx2_arr
