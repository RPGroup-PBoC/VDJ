"""
Parameter Estimation for PC Leaving Rate
--------------------------------------------------------------------------------
Author: Griffin Chure
Last Modified: September 25, 2019
License: MIT
Associated Files: `expon_model.stan`

Description
--------------------------------------------------------------------------------
This script performs parameter inference of the leaving rate (parameterized by tau) given the dwell time distributions. 

Notes
--------------------------------------------------------------------------------
This script is designed to be exectued in the `code/analysis` directory and
loads the necessary data via a relative path.
"""
import numpy as np
import pandas as pd
import tqdm
import vdj.io
import vdj.bayes
import vdj.stats

# Load data and stan model
data = pd.read_csv('../../data/compiled_dwell_times.csv', comment='#')
model = vdj.bayes.StanModel('../stan/expon_dwell_model.stan', 
                            force_compile=True)

# Iterate through the data and fit while storing thinned samples
samps_df = []
stats_df = []
for g, d in tqdm.tqdm(data.groupby(['mutant', 'salt'])):
    _, samples = model.sample({'N':len(d), 'dwell':d['dwell_time_min']},
                            control=dict(adapt_delta=0.95), iter=5000)
    stats = model.summary()

    # Parse the mutant
    mut = vdj.io.mutation_parser(g[0])

    # Add identifiers and append
    samples['mutant'] = g[0]
    samples['seq'] = mut['seq']
    samples['n_muts'] = mut['n_muts']
    samples['salt'] = g[1]
    stats['mutant'] = g[0]
    stats['seq'] = mut['seq']
    stats['n_muts'] = mut['n_muts']
    stats['salt'] = g[1]
    samps_df.append(samples.iloc[::10])
    stats_df.append(stats)

# Concatenate and save the dataframes
samps = pd.concat(samps_df)
samps.to_csv('../../data/exponential_fitting_samples.csv', index=False)
stats = pd.concat(stats_df)
stats.to_csv('../../data/exponential_fitting_summary.csv', index=False)