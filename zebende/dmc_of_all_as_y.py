import numpy as np
from . import ordering_x_dmcx2_of



def dmc_of_all_as_y(data):
    aux = list(range(data.shape[1]))
    dmc_list = [aux]
    for i in range(data.shape[1]-1):
        aux = aux[1:] + aux[:1]

        

        
        dmc_list.append(aux)

    dmc_list = np.array(dmc_list)

    dmc_list = ordering_x_dmcx2_of(dmc_list)

    return dmc_list

