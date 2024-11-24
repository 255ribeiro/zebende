# %%
import numpy as np
import zebende as zb
import pandas as pd

# %%
tws = pd.read_csv('./test_data/time_window_scales.txt', header=None).to_numpy().flatten()
tws

# %%
data_1 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_F3.txt', header=None)
data_2 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_F6.txt', header=None)
data_3 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_P3.txt', header=None)
data_4 = pd.read_csv('./test_data/data_pre_proc/S001/S001_03/S001_03_P6.txt', header=None)

# %%
data = pd.concat([data_1, data_2, data_3, data_4], axis=1)
data.to_csv('./test_data/data_pre_proc/S001/four_channels.csv')

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
