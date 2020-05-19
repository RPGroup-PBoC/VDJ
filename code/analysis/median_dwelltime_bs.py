#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the dwell time data
dwell = pd.read_csv('../../data/compiled_dwell_times.csv')

# Ignore the other HMGB1 concentrations
dwell = dwell[dwell['hmgb1']==80].copy()

#%%
# Define a bootstrapping function
def bs_median(dist, rep=1E6):
    n_meas = len(dist)
    draws = np.random.choice(dist, size=(n_meas, int(rep)), replace=True)
    medians = np.mean(draws, axis=0)
    # percs = np.percentile(([1, 2.5, 25, 45, 55, 75, 97.5, 99]))
    # _df = pd.DataFrame([])
    return medians

#%%
# Look at one sample
samp = dwell[dwell['mutant']=='12HeptA6T']
bs = bs_median(samp['dwell_time_min'])

#%%
bs_median(samp['dwell_time_min'].values)

#%%
