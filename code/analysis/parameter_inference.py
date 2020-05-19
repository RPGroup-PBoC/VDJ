# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import vdj.bayes
import vdj.stats
import tqdm

# Load the data sets
dwell = pd.read_csv('../../data/compiled_dwell_times.csv')
fates = pd.read_csv('../../data/compiled_cutting_events.csv')
f_looped = pd.read_csv('../../data/compiled_looping_fraction.csv')


# Begin with the pooled analysis to get a feeling for the state of the data
model = vdj.bayes.StanModel('../stan/pooled_model.stan') 

# Instantiate empty lists to add dataframes
samples_dfs, stats_dfs = [], []
for i, m in enumerate(tqdm.tqdm(f_looped['mutant'].unique(), 
                       desc='Pooled model inference')):
    # Get the data for each mutant.
    _dwell = dwell[dwell['mutant']==m]
    _fates = fates[fates['mutant']==m]
    _f_looped = f_looped[f_looped['mutant']==m]

    # Instantiate the data dictionary. 
    data_dict = {'N':len(_dwell), 'n_beads':int(_fates['n_beads'].sum()), 
                 'n_cuts':int(_fates['n_cuts'].sum()), 
                 'total_frames':int(_f_looped['total_frames'].sum()),
                 'looped_frames':int(_f_looped['looped_frames'].sum()),
                 'dwell_time':_dwell['dwell_time_min']}

    # Sample the model and compute the statistics
    _, samples = model.sample(data_dict)
    stats = model.summary()
    samples['mutant'] = m
    stats['mutant'] = m
    samples_dfs.append(samples)
    stats_dfs.append(stats)

# Concatenate the data frames and save to disk
samples = pd.concat(samples_dfs)
stats = pd.concat(stats_dfs)
samples.to_csv('../../data/pooled_model_samples.csv', index=False)
stats.to_csv('../../data/pooled_model_summary.csv', index=False)

