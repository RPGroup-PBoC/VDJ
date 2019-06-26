# -*- coding: utf-8 -*-
#%% 
import numpy as np
import pandas as pd
from scipy.special import gammaln, logsumexp
import vdj.io

# Define the posterior distribution
def log_posterior(N, n, p):
    binom_coeff = gammaln(N + 2) - gammaln(n + 1) - gammaln(N - n + 1)
    bernoulli = n * np.log(p) + (N - n) * np.log(1 - p)
    return binom_coeff  + bernoulli

# Load the fates
data = pd.read_csv('../../data/compiled_bead_fates.csv')

# Get the mutant info
mut_info = {m:vdj.io.mutation_parser(m) for m in data['mutant'].unique()}

# Compute the pooled statistics
pooled = data.groupby(['mutant', 'date', 'salt', 'hmgb1']).agg(('sum')).reset_index()
pooled = pooled[['mutant', 'salt', 'date', 'hmgb1', 'n_beads', 'n_cuts']]
pooled['mode'] = pooled['n_cuts'].values / pooled['n_beads']
pooled['std'] = (pooled['n_cuts'].values * (pooled['n_beads'] -\
                 pooled['n_cuts'])) / pooled['n_beads'].values**3
for m, seq in mut_info.items():
    pooled.loc[pooled['mutant']==m, 'n_muts'] = seq['n_muts']
pooled.to_csv('../../data/pooled_cutting_probability.csv', index=False)

# Compute the hierarchical stats
hier = data.copy()
hier['mode'] = hier['n_cuts'].values / hier['n_beads']
hier['std'] = (hier['n_cuts'].values * (hier['n_beads'] -\
             hier['n_cuts'].values)) / hier['n_beads']**3
for m, seq in mut_info.items():
    hier.loc[hier['mutant']==m, 'n_muts'] = seq['n_muts']
hier.to_csv('../../data/independent_cutting_probability.csv', index=False)

#%%
# Compute the posteriors for the pooled statistics
prob_range = np.linspace(1E-5, 1 - 1E-5, 200)
posts = []

for g, d in pooled.groupby(['mutant', 'salt', 'hmgb1']):
    # Evaluate the posterior over the full range of probabilities
    N = d['n_beads'].values[0]
    n = d['n_cuts'].values[0]
    log_post = log_posterior(N, n, prob_range) 
    _post = log_post - logsumexp(log_post) 
    _post = np.exp(_post)

    # Store the results in a dataframe
    _df = pd.DataFrame(np.array([prob_range, _post]).T, columns=['probability', 
                                                                  'posterior'])
    _df['mutant'] = g[0]
    _df['salt'] = g[1]
    _df['hmgb1'] = g[2]
    _df['n_muts'] = d['n_muts'].unique()[0]
    posts.append(_df)

# Concatenate the posteriors and save to disk
pooled_posts = pd.concat(posts).reset_index()
pooled_posts.to_csv('../../data/pooled_cutting_probability_posteriors.csv', 
                   index=False)



_post

#%%


#%%
