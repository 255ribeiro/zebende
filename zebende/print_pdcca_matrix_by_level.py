import numpy as np
from numpy.typing import NDArray


def print_pdcca_matrix_by_level(pdcca:NDArray[np.float64]):
    """Print the <span>&Rho;<sub>DCCA</sub></span> matrix by level.

    Args:
        pdcca (NDArray[np.float64]): <span>&Rho;<sub>DCCA</sub></span> calculation in the matrix format.
    """

    print(pdcca.transpose(2,0,1))