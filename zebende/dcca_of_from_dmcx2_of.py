from itertools import combinations
import numpy as np


def dcca_of_from_dmcx2_of(dmc_of):
    if type(dmc_of) == np.ndarray:
        dmc_of = dmc_of.tolist()
    out = []
    for i in range(len(dmc_of)):
        temp = list(combinations(dmc_of[i],2))
        for j in temp:
            j = list(j)

            j.sort()

            if not(j in out):
                out.append(j)
    return np.array(out)