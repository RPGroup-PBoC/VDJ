#%%
import numpy as np
import pandas as pd 
import numba
from tqdm import tqdm


# Define the number of bootstrap replicates
n_bs = int(1E7)

# Load the cutting data 
data = pd.read_csv('../../data/compiled_bead_fates.csv')
data = data[(data['hmgb1']==80) & (data['salt']=='Mg')]

# Isolate the wt12rss data and compute quantities
ref_data = data[(data['mutant']=='WT12rss')]
ref_beads = np.zeros(int(ref_data['n_beads'].sum()))
ref_beads[:int(ref_data['n_cuts'].sum())] = 1
ref_pcut = ref_beads.sum() / len(ref_beads)


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
        exp = concat[:len(x)].sum() / len(x)
        ctrl = concat[len(x):].sum() / len(y) 
        out[i] = np.abs(exp- ctrl)
    return out
#%%
p_val = {}
for g, d in data[data['mutant']!='WT12rss'].groupby(['mutant']):
    print(f'Processing mutant {g}...')
    # Get the number of loops  and total number of observations
    if d['n_beads'].sum() != 0:
        exp_beads = np.zeros(int(d['n_beads'].sum()))
        exp_beads[:int(d['n_cuts'].sum())] = 1
        exp_pcut = exp_beads.sum() / len(exp_beads)

        # Compute the observed difference in cut probability
        diff_pcut = np.abs(exp_pcut - ref_pcut)

        # Perform the simulation
        out = permute(exp_beads, ref_beads)

        # Compute the p-value
        p_val[g] = np.sum(out >= diff_pcut) / len(out)

# %%
keys = list(p_val.keys())
vals = list(p_val.values())
df = pd.DataFrame(np.array([keys, vals]).T, columns=['mutant', 'p_value'])
df.to_csv('../../data/cutting_probability_p_values.csv', index=False)
