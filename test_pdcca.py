# %%
import numpy as np
import zebende as zb
import pandas as pd

# %%
tws = pd.read_csv('./test_data/time_window_scales.txt', header=None).to_numpy().flatten()
tws =  np.array([4, 5, 7, 9, 11, 13, 16, 20, 23, 28, 33, 38, 45, 52, 60, 69, 79, 91, 104, 119, 135, 154, 174, 198, 223, 252, 285, 321,
  362, 407, 457, 513, 575, 645, 723, 809, 905, 1011, 1130, 1261, 1407, 1570])
print(tws)
# %%
data_1 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_F3.txt', header=None, names=['F3'])
data_2 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_F6.txt', header=None, names=['F6'])
data_3 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_P3.txt', header=None, names=['P3'])
data_4 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_P6.txt', header=None, names=['P6'])

# %%
data = pd.concat([data_1, data_2, data_3, data_4], axis=1)
data.to_csv('./test_data/data_pre_proc/S001/S001_T03_F3_F6_P3_P6.csv', index=False)

# %%
data.tail()

# %%
data = data.to_numpy()
data

# %%
int_data = zb.integrated_series(data)
int_data.shape[0]


# %%
dfa , dcca, pdcca = zb.p_dcca(int_data, tws)
