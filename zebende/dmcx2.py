from typing import Literal

import numpy as np

from . import (
    dcca_of_from_dmcx2_of,
    dmc_of_all_as_y,
    ordering_x_dmcx2_of,
    p_dcca,
)

ENUM_DMCx2_of = Literal['all-full', 'first-full']

def dmcx2(data: np.ndarray, tws: np.ndarray, dmcx2_of: np.ndarray | list | ENUM_DMCx2_of = 'all-full', time_steps: np.ndarray | None = None, DCCA_of: np.ndarray | list | None = None):
    
    if type(dmcx2_of) == str:
        # creating ndarray of y and x values for DMCx2 calculations
        if dmcx2_of == 'first-full':
            dmcx2_of =np.array( [np.arange(data.shape[1])])
        # creating ndarray of y and x values for DMCx2 calculations
        elif dmcx2_of == 'all-full':
            dmcx2_of = dmc_of_all_as_y(data)
    else:
        dmcx2_of = ordering_x_dmcx2_of(dmcx2_of)

    # creating ndarray for P_DCCA calculations based on the DMCx2 array
    if DCCA_of == None:
        DCCA_of = dcca_of_from_dmcx2_of(dmcx2_of)
    # P_DCCA calculations
    F_DFA_arr, DCCA_arr, P_DCCA_arr = p_dcca(data=data, tws=tws, time_steps=time_steps, DCCA_of=DCCA_of, P_DCCA_output_format = 'matrix')

    # DMCx2 output matrix
    DMCx2_arr = np.empty(shape=(tws.shape[0], dmcx2_of.shape[1]), dtype=np.float64)

    for n_index in range(len(tws)):

        for j in range(dmcx2_of.shape[0]):

            y_indexes = dmcx2_of[j, 0:1]
            x_indexes = dmcx2_of[j, 1:]

            mat_x = P_DCCA_arr[np.ix_(x_indexes, x_indexes)][:, :, n_index]
            vec_y = P_DCCA_arr[np.ix_(y_indexes, x_indexes)][:, :, n_index]

            # dmcx2 calculation
            DMCx2_arr[n_index, j] = vec_y @ np.linalg.inv(mat_x) @ vec_y.T

    return F_DFA_arr, DCCA_arr, P_DCCA_arr, DMCx2_arr