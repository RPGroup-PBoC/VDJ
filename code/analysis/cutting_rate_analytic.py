# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd 
import vdj.io
import vdj.model
import vdj.stats
import imp
imp.reload(vdj.model)
imp.reload(vdj.stats)

# Load the data 
data = pd.read_csv('../../data/compiled_dwell_times.csv')
# %%
# Instantiate the dataframe 
percs_dfs = pd.DataFrame([])

# Define the range of constants overwhich to query 
tau = np.linspace(1, 500, 5000)

# Loop through each mutant and evaluate the log posterior
stat_df = pd.DataFrame()
post_dfs = []
for g, d in data.groupby('mutant'):
    cuts = d[d['cut']==1]['dwell_time_sec'].values / 60
    unloops = d[d['cut']==0]['dwell_time_sec'].values / 60
    log_post = vdj.model.cutting_rate_log_posterior(tau, cuts, unloops)

    # Rescale the log posterior
    log_post -= log_post.max() 

    # Compute the cdf and percentiles
    cdf = np.cumsum(np.exp(log_post))
    cdf /= cdf.max()
    percs = vdj.stats.compute_percentiles(cdf, tau)

    # Extract the interesting variables
    low, high = percs[0.95]
    width = high - low
    median = percs['median']

    # Get the mutant identifying information. 
    seq = vdj.io.mutation_parser(g)
    stat_df = stat_df.append({'median':median, '95_low':low, '95_high':high, 'width':width, 
                    'mutant':g, 'seq':seq['seq'], 'n_muts':seq['n_muts']}, ignore_index=True)

    _post = pd.DataFrame(np.array([log_post, 
        np.exp(log_post) / np.sum(np.exp(log_post)), 
        cdf, tau]).T, columns=['log_posterior', 'posterior_pdf', 'posterior_cdf', 'tau'])
    _post['mutant'] = g
    _post['seq'] = seq['seq']
    _post['n_muts'] = seq['n_muts']
    post_dfs.append(_post)
post_df = pd.concat(post_dfs)

# %%
# Save things to disk
post_df.to_csv('../../data/cutting_rate_posteriors.csv', index=False)
stat_df.to_csv('../../data/cutting_rate_analytic_summary.csv', index=False)


