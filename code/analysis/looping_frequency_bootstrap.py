#%%
import numpy as np
import pandas as pd
import vdj.io
import vdj.stats
from tqdm import tqdm

#%%
data = pd.read_csv('../../data/compiled_looping_events.csv', sep=',')

#%%
# Need to group by mutant, Ca and HMGB1 concentration
bs_df = pd.DataFrame(columns=['mean','bs_low','bs_high',
                            'mutant','salt','hmgb1'])
mutant_count=0
for g,d in data.groupby(['mutant', 'salt', 'hmgb1']):
    mutant_count += 1
    print(mutant_count)
    sampling = np.zeros(1000000)
    for n in range(0, 1000000, 1):
        sampling[n] = np.sum(d['n_loops'].sample(n=len(d),replace=True)) / len(d)
    bs_low, bs_high = vdj.stats.compute_hpd(sampling, mass_frac=0.95)
    mean = np.mean(sampling)
    df = pd.DataFrame({'mean':[mean], 
                        'bs_low':[bs_low], 
                        'bs_high':[bs_high]})
    df['mutant'] = g[0]
    df['salt'] = g[1]
    df['hmgb1'] = g[2]
    bs_df = bs_df.append(df)

#%%
