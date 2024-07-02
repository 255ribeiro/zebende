import numpy as np


def dmc_of_all_as_y(data):
    temp = list(range(data.shape[1]))
    dmc_list = []
    for i in range(data.shape[1]):
        temp = temp[i:] + temp[:i]
        dmc_list.append(temp)

if __name__ == '__main__':
    ...
    data = np.arange(16*16).reshape(16,16)
    print(data)
    test_01 = dmc_of_all_as_y(data)
    test_02 = np.arange(data.shape[1])
    print(test_01)
    print(test_02)