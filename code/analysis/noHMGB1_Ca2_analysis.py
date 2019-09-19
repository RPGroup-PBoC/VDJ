#%%[markdown]
# # Experiments to back out looping and cutting rates

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.io
import vdj.viz
import vdj.stats
vdj.viz.plotting_style()
imp.reload(vdj.io)

#%% [markdown]
# Soichi did a few experiments different than normal to see if we could back out
# some of the parameters in isolation. One of these experiments was measuring
# the looping in the absence of HMGB1 while another was doing the whole thing
# with Ca2+ instead of Mg2+ which should inhibit cutting. A concern with the
# latter experiment is that subbing out Ca inplace of Mg may alter the
# flexibility of the DNA and therefore change the looping frequency. This may be
# worth talking to Rob and or SLJ about as I'm sure they played with this
# before. FOr now, we will just approach the first case in which HMGB1 was
# absent. We'll start with loading the mat file and extracting the data. 

#%%
mat_ca = vdj.io.ProcessTPM('../../data/12SpacG11T_Ca_analysis_270_analyzed.mat')
mat_mg = vdj.io.ProcessTPM('../../data/12SpacG11T_analysis_280_analyzed.mat')

#%%
fit, sample_ca, stats_ca = mat_ca.run_inference(stan_model='../stan/hierarchical_model.stan')
fit, sample_mg, stats_mg = mat_mg.run_inference(stan_model='../stan/hierarchical_model.stan')

#%%
dwell_ca = mat_ca.dwell
dwell_mg = mat_mg.dwell

#%%
r_ca = vdj.stats.compute_hpd(sample_ca['k_unloop'].values, 0.95)
r_mg = vdj.stats.compute_hpd(sample_mg['r_cut'] + sample_mg['k_unloop'].values, 0.95)
f_range = np.logspace(0, 5)
low = 1 - np.exp(-r_ca[0] * f_range)
high = 1 - np.exp(-r_ca[1] * f_range)

# %%
plt.hist(dwell['dwell_time_ms'], bins=75, density=True)
plt.plot(f_range, r * np.exp(-r * f_range), '-', color='tomato')
plt.plot(f_range, r * np.exp(-r * f_range), '-', color='tomato')

#%%
# plt.fill_between(f_range, low, high, color='tomato', alpha=0.5)
plt.step(np.sort(dwell_mg['dwell_time_ms'].values), np.linspace(0, 1, len(dwell_mg)),
color='tomato')
plt.step(np.sort(dwell_ca['dwell_time_ms'].values), np.linspace(0, 1, len(dwell_ca)),
color='dodgerblue')


#%%


#%%
