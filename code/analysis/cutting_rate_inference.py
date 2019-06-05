# -*-  coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import tqdm
import vdj.io
import vdj.bayes

# %%
# Load teh dwell time data sets.
dwell = pd.read_csv('../../data/compiled_dwell_times.csv')

# Load the stan model for the pooled case
model = vdj.bayes.StanModel('../stan/pooled_cutting_rate.stan', 
                            force_compile=True)

# %% Perform the inference
# Set up the storage list for the summaries
stats = []

# Iterate through each mutant. 
for g, d in tqdm.tqdm(dwell.groupby('mutant')):
    # Define the data dictionary
    cuts = d[d['cut']==1]
    unloops = d[d['cut']==0]
    data_dict = {'N':len(cuts), 'M':len(unloops), 
                'cut':cuts['dwell_time_s'], 'unloop':unloops['dwell_time_s']}
    # Sample and compute the parameter summary
    model.sample(data_dict)
    summary = model.summary(parnames=['tau', 'r'])
    summary['mutant'] = g

    # Get the sequences
    seq = vdj.io.mutation_parser(g)
    summary['seq'] = seq['seq']
    summary['n_muts'] = seq['n_muts']
    stats.append(summary)
stats = pd.concat(stats)

# %% 
# Save it to CSV 
stats.to_csv('../../data/pooled_cutting_rate.csv', index=False)

#%%
dwell

#%%
stats

#%%
