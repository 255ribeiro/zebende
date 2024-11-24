import numpy as np

def print_pdcca_matrix_by_level(pdcca):
    print(pdcca.transpose(2,0,1))