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
    sampling = np.random.choice(d['n_loops'].values,size=(len(d),10**6),replace=True)
    loop_freq = np.sum(sampling, axis=1) / len(d)
    bs_low, bs_high = vdj.stats.compute_hpd(loop_freq, mass_frac=0.95)
    mean = np.mean(sampling)
    df = pd.DataFrame({'mean':[mean], 
                        'bs_low':[bs_low], 
                        'bs_high':[bs_high]})
    df['mutant'] = g[0]
    df['salt'] = g[1]
    df['hmgb1'] = g[2]
    bs_df = bs_df.append(df, ignore_index=True)
bs_df = bs_df[['mutant', 'salt', 'hmgb1', 'mean', 'bs_low', 'bs_high']]
#%%
bs_df.to_csv('../../data/compiled_loop_freq_bs.csv')

#%%
