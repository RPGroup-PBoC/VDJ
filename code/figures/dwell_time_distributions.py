#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import vdj.io 
import vdj.viz
vdj.viz.plotting_style()

#%%
# load the dwell time data
data = pd.read_csv('../../data/compiled_dwell_times.csv')
data = data[(data['salt']=='Mg') & (data['hmgb1']==80)]
ref = d[d['mutant']=='WT12rss']

# Load the fitting statistics. 
samps = pd.read_csv('../../data/expon_samples.csv')
samps = samps[(samps['mutant']=='WT12rss') & (samps['salt']=='Mg')]
stats = pd.read_csv('../../data/expon_summary.csv')
stats = stats[(stats['mutant']=='WT12rss') & (stats['salt']=='Mg')]

#%%
# Set up the figure canvas
fig, ax = plt.subplots(1, 2, figsize=(4.5, 2))
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
ax[1].set_xscale('log')
ax[0].set_yscale('log')

ax[0].hist(ref['dwell_time_min'], density=True, bins=30, color='slategrey', alpha=0.75)
x = np.sort(ref['dwell_time_min'])
y = np.arange(0, len(x), 1) / len(x)
ax[1].step(x, y, '-', color='w', lw=2, label='__nolegend__')
ax[1].step(x, y, '-', color='slategrey', lw=1.25, )
#%%
