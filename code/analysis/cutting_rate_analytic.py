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
tau = np.linspace(0, 100, 500)

# Loop through each mutant and evaluate the log posterior
for g, d in data.groupby('mutant'):
    cuts = d[d['cut']==1]['dwell_time_min'].values
    unloops = d[d['cut']==0]['dwell_time_min'].values
    log_post = vdj.model.cutting_rate_log_posterior(tau, cuts, unloops)

    # Compute the cdf and percentiles
    cdf = np.cumsum(np.exp(log_post))
    cdf /= cdf.sum()
    percs = vdj.stats.compute_percentiles(cdf, tau)
    break


# %%
vdj.model.cutting_rate_log_posterior()


#%%
log_post


#%%
log_post

#%%

# %%
def cutting_rate_log_prior(tau, alpha=0.86, beta=9.03):
    return st.invgamma(alpha, scale=beta, loc=0).logpdf(tau).sum()

def cutting_rate_log_posterior(tau, t_cut, t_unloop, **kwargs):
    log_like = cutting_rate_log_likelihood(tau, t_cut, t_unloop)
    log_prior = cutting_rate_log_prior(tau, **kwargs)
    return log_like + log_prior



