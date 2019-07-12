# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import vdj.bayes
import vdj.stats
import tqdm 

# Load the dwell time data
data = pd.read_csv('../../data/compiled_dwell_times.csv')

# Only consider HMGB1 at 80nM
data = data[data['hmgb1']==80]

# Identify mutants for which calcium and magnesium data exist.
muts = [g for g, d in data.groupby('mutant') if len(d['salt'].unique()) == 2]

# Isolate only those mutants in the data set. 
#%%
salt_data = data[(data['mutant']==muts[0]) | (data['mutant']==muts[1]) |
                 (data['mutant']==muts[2])]
#%%
# Load the stan model used for fitting.
model = vdj.bayes.StanModel('../stan/expon_dwell_model.stan')

#%% Set up the storage sample and posterior predictive data frames. 
samples_dfs = []
ppc_dfs = []

# Define the deadfilter for cdf correction
dead_filter = 21 / 60 # 21 seconds / 60 seconds => min
# Iterate through each salt and mutant pairing. 
for g, d in tqdm.tqdm(salt_data.groupby(['mutant', 'salt']), desc='Mutant and salts'):
    # Define the data dictionary. 
    data_dict = {'N': len(d), 'dwell':d['dwell_time_min']}

    # Sample the model. 
    _, samples = model.sample(data_dict)
    samples['mutant'] = g[0]
    samples['salt'] = g[1]
    samples_dfs.append(samples)

    # Perform posterior predictive checks
    _dfs = []
    for i, s in enumerate(samples['tau'].values):
        draws = np.random.exponential(s, size=data_dict['N'])
        # Make sure there aren't any below the dead time
        n_low = np.sum(draws < dead_filter)
        if n_low > 0:
            draws = list(draws[draws > dead_filter])
            while n_low > 0:
                _draw = np.random.exponential(s, size=int(n_low))
                for _d in _draw[_draw > dead_filter]:
                    draws.append(_d)
                n_low = np.sum(_draw < dead_filter)
        _df = pd.DataFrame([])
        _df['draws'] = np.sort(draws)
        _df['ecdf'] = np.arange(0, data_dict['N']) / data_dict['N']
        _df['sample_idx'] = i + 1
        _dfs.append(_df)

    # Assemble the posterior predictive checks into a single data frame
    _df = pd.concat(_dfs)
    _df['mutant'] = g[0]
    _df['salt'] = g[1]
    ppc_dfs.append(_df)

# Generate a final posterior predictive data frame and save to disk.
ppc = pd.concat(ppc_dfs)
ppc.to_csv('../../data/CaMg_dwelltime_exponential_ppc.csv', index=False)


#%%
