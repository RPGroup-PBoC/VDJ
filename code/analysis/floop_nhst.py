#%%
import numpy as np
import pandas as pd 
import numba
from tqdm import tqdm


# Define the number of bootstrap replicates
n_bs = int(1E7)

# Load the cutting data 
data = pd.read_csv('../../data/compiled_looping_events.csv')
data = data[(data['hmgb1']==80) & (data['salt']=='Mg')]

# Islate the wt12rss data. 
ref_data = data[(data['mutant']=='WT12rss')]
ref_loops = ref_data['n_loops'].values

# Compute the looping frequency of the reference sequence. 
ref_loop_freq = ref_data['n_loops'].sum() / len(ref_data)

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
        freq = concat[:len(x)].sum() / len(x)
        ctrl = concat[len(x):].sum() / len(y) 
        out[i] = np.abs(freq - ctrl)
    return out

# %% Group by mutant
p_val = {}
for g, d in data[data['mutant']!='WT12rss'].groupby(['mutant']):
    print(f'Processing mutant {g}...')
    # Get the number of loops  and total number of observations
    n_loops = d['n_loops'].values
    n_obs = len(d)

    # Compute the observed difference in loop frequency
    loop_freq_diff = np.abs((n_loops.sum() / n_obs) - ref_loop_freq)

    out = permute(n_loops, ref_loops)

    # Compute the p-value
    p_val[g] = np.sum(out >= loop_freq_diff) / len(out)

# %%
keys = list(p_val.keys())
vals = list(p_val.values())
df = pd.DataFrame(np.array([keys, vals]).T, columns=['mutant', 'p_value'])
df.to_csv('../../data/looping_frequency_p_values.csv', index=False)

# %%
df


# %%
