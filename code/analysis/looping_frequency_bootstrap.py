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
col_names = [("bs_95_low", "bs_95_high"), ("bs_75_low", "bs_75_high"),
             ("bs_50_low", "bs_50_high"),  ("bs_25_low", "bs_25_high"),
             ("bs_10_low", "bs_10_high"), ("bs_5_low", "bs_5_high")]
bs_reps = int(1E6)
bs_df = pd.DataFrame([])
for g,d in tqdm.tqdm(data.groupby(['mutant', 'salt', 'hmgb1'])):
    sampling = np.random.choice(d['n_loops'].values,size=(len(d), bs_reps),replace=True)
    loop_freq = np.sum(sampling, axis=0) / len(d)
    df_dict = {'mutant':g[0], 'salt':g[1], 'hmgb1':g[2],
                'n_loops':d['n_loops'].sum(), 'n_beads':len(d),
                'loops_per_bead':d['n_loops'].sum() / len(d)}
    for perc, col in zip(percentiles, col_names):
        computed_percentiles = np.percentile(loop_freq, perc)
        for i in range(2):
            df_dict[col[i]] = computed_percentiles[i]

    bs_df = bs_df.append(df_dict, ignore_index=True)
bs_df.to_csv('../../data/compiled_loop_freq_bs.csv')

#%%
