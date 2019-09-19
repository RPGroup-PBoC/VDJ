#%%
import numpy as np
import pandas as pd
import vdj.io
import vdj.stats
import tqdm

#%%
data = pd.read_csv('../../data/compiled_looping_events.csv')

#%%
percentiles = [(2.5, 97.5), (87.5, 12.5), (25, 75), (37.5, 62.5), (45, 55), 
               (47.5, 52.5)]
percentile_names = [95, 75, 50, 25, 10, 5]

bs_reps = int(1E6)
bs_df = pd.DataFrame([])
for g,d in tqdm.tqdm(data.groupby(['mutant', 'salt', 'hmgb1'])):
    sampling = np.random.choice(d['n_loops'].values,size=(len(d), 
                                bs_reps),replace=True)
    loop_freq = np.sum(sampling, axis=0) / len(d)
    for perc, name in zip(percentiles, percentile_names):
        computed = np.percentile(loop_freq, perc)
        df_dict = {'mutant':g[0], 'salt':g[1], 'hmgb1':g[2],
                'n_loops':d['n_loops'].sum(), 'n_beads':len(d),
                'loops_per_bead':d['n_loops'].sum() / len(d),
                'percentile': name, 'low':computed[0], 
                'high':computed[1]} 
        bs_df = bs_df.append(df_dict, ignore_index=True)

bs_df.to_csv('../../data/compiled_loop_freq_bs.csv')

#%%
