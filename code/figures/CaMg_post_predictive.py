# -*- coding: utf-8 -*- 
#%%
import numpy as numpy
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import vdj.viz
vdj.viz.plotting_style()

# Load the posterior predictive samples
ppc = pd.read_csv('../../data/CaMg_dwelltime_exponential_ppc.csv')

# Load the actual data
data = pd.read_csv('../../data/compiled_dwell_times.csv')

# %%
# Set up the figure canvas
fig = plt.figure(figsize=(7.5, 4))
gs = GridSpec(2, 6)

# Assign the axes
ppc_ax = [fig.add_subplot(gs[0, i]) for i in range(6)]
ax0 = fig.add_subplot(gs[1, :2])
ax1 = fig.add_subplot(gs[1, 2:4])
ax2 = fig.add_subplot(gs[1, 4:])

# Set the xscaling to log
for a in ppc_ax:
    a.set_xscale('log')

# Format the labels
for i, a in enumerate(ppc_ax):
    if i > 0:
        a.set_yticks([])
    if i%2 == 0:
        a.set_title('Ca$^{2+}$')    
    else:
        a.set_title('Mg$^{2+}$')

# Define the mutant axes:
mut_ax = {'WT12rss':ppc_ax[0:2], '12HeptA4T':ppc_ax[2: 4], '12SpacG11T':ppc_ax[4:]}
salt_ax = {'Ca': 0, 'Mg': 1}
colors = {'Ca': 'green', 'Mg':'purple'}

# Iterate through the mutants and salts and plot the data. 
for m, a in mut_ax.items():
    for s, ad in salt_ax.items():
        d = data[(data['mutant']==m) & (data['salt']==s)]
        x = np.sort(d['dwell_time_min'].values)
        y = np.arange(0, len(x)) / len(x)
        a[ad].step(x, y, 'k-', lw=0.75)

        # Plot the PPC
        _ppc = ppc[(ppc['mutant']==m) & (ppc['salt']==s)]
        for g, d in _ppc.groupby(['sample_idx']):
            if g%100 == 0:
                a[ad].step(d['draws'] - 21/60, d['ecdf'], '-', lw=0.05, alpha=0.5, color=colors[s])

#%%
