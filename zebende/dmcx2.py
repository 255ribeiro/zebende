import numpy as np

from . import p_dcca
from . import dcca_of_from_dmcx2_of
from . import dmc_of_all_as_y
from . import ordering_x_dmcx2_of


def dmcx2(data: np.ndarray, tws: np.ndarray, dmcx2_of: np.ndarray| list | str = 'all-full', time_steps: np.ndarray | None = None, DCCA_of: np.ndarray | list | None = None):
    # creating ndarray of y and x values for DMCx2 calculations 
    if dmcx2_of == 'first-full':
        dmcx2_of = np.arange(data.shape[1])
    # creating ndarray of y and x values for DMCx2 calculations 
    elif dmcx2_of=='all-full':
        dmcx2_of = dmc_of_all_as_y(data)
    else:
        dmcx2_of = ordering_x_dmcx2_of(dmcx2_of)
    # creating ndarray for P_DCCA calculations based on the DMCx2 array
    if DCCA_of == None:
        DCCA_of = dcca_of_from_dmcx2_of(dmcx2_of)
    # P_DCCA calculations             
    F_DFA_arr, DCCA_arr, P_DCCA_arr = p_dcca(data=data, tws=tws, time_steps=time_steps, DCCA_of=DCCA_of)

    # DMCx2 output matrix
    DMCx2_arr = np.empty(shape=(tws.shape[0], dmcx2_of.shape[1]), dtype=np.float64)
    # DMcx2 auxilary matrix
    dmcx2_mat = np.ones(shape=(dmcx2_of.max(), dmcx2_of.max(), tws.shape[0]), dtype=np.float64)

    for n_index in range(len(tws)):
        # populating the matix
        for i in range(len(DCCA_of.shape[1])):
            dmcx2_mat[DCCA_of[i, 0],DCCA_of[i, 1], n_index] = P_DCCA_arr[n_index, i]
            dmcx2_mat[DCCA_of[i, 1],DCCA_of[i, 0], n_index] = P_DCCA_arr[n_index, i]
            # extracting values for calculations
            for j in range(dmcx2_of.shape[1]):
                y_indexes = dmcx2_of[j, 0]
                x_index = dmcx2_of[j, 1:]
                mat_x = dmcx2_mat[y_indexes, y_indexes, n_index]
                vec_y = dmcx2_mat[x_index, y_indexes, n_index]
                # dmcx2 calculation
                DMCx2_arr[n_index, j] = vec_y @ np.linalg.inv(mat_x) @ vec_y

    return F_DFA_arr, DCCA_arr, P_DCCA_arr, DMCx2_arr
