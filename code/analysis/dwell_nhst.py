#%%
import numpy as np 
import pandas as pd 

# Define the number of permutations. 
n_bs = int(1E6)

# Load the data and restrict to useful conditions
data = pd.read_csv('../../data/compiled_dwell_times.csv')
data = data[(data['salt']=='Mg') & (data['hmgb1']==80)]

# Find the reference distribution and compute the median dwell time. 
ref_dwell = data[data['mutant']=='WT12rss']['dwell_time_min'].values
median_dwell = np.median(ref_dwell)

#%%
# Define functions for the NHST with numba for speed
@numba.njit 
def draw(x, y):
    concat = np.concatenate((x,y))
    np.random.shuffle(concat)
    return concat

@numba.njit
def permute(x, y, size=n_bs):
    out = np.empty(size)
    for i in range(size):
        concat = draw(x, y) 
        exp_median = np.median(concat[:len(x)])
        ctrl_median = np.median(concat[len(x):])
        out[i] = np.abs(exp_median - ctrl_median)
    return out


# %% Group by mutant
p_val = {}
for g, d in data[data['mutant']!='WT12rss'].groupby(['mutant']):
    print(f'Processing mutant {g}...')
    # Get the number of loops  and total number of observations
    dwell = d['dwell_time_min'].values

    # Compute the observed difference in loop frequency
    dwell_diff = np.abs(np.median(dwell) - median_dwell)

    out = permute(dwell, ref_dwell)

    # Compute the p-value
    p_val[g] = np.sum(out >= dwell_diff) / len(out)

# %%
keys = list(p_val.keys())
vals = list(p_val.values())
df = pd.DataFrame(np.array([keys, vals]).T, columns=['mutant', 'p_value'])
df.to_csv('../../data/dwell_time_p_values.csv', index=False)


# %%
