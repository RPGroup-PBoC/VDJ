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
fig = plt.figure(figsize=(7.5, 4.3))
gs = GridSpec(2, 6)

# Assign the axes
ppc_ax = [fig.add_subplot(gs[0, i]) for i in range(6)]
ax0 = fig.add_subplot(gs[1, :2])
ax1 = fig.add_subplot(gs[1, 2:4])
ax2 = fig.add_subplot(gs[1, 4:])

# Set the xscaling to log
for a in ppc_ax:
    a.set_xscale('log')
    a.set_xticks([100, 1000])
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
for i, a in enumerate([ax0, ax1, ax2]):
    a.set_xscale('log')
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xticks([10, 100, 1000, 10000])
    a.set_xlabel('dwell time [s]', fontsize=8)
    a.set_title('dwell time [s]', fontsize=8, y=1.04)
    if i > 0:
        a.set_yticks([])

# Format the labels
for i, a in enumerate(ppc_ax):
    if i > 0:
        a.set_yticks([])
    if i%2 == 0:
        a.set_title('Ca$^{2+}$', fontsize=8)   
    else:
        a.set_title('Mg$^{2+}$', fontsize=8)

# Add axis y axis label
ppc_ax[0].set_ylabel('cumulative distribution', fontsize=8)
ax0.set_ylabel('cumulative distribution', fontsize=8)

# Set the panel labels
fig.text(0.21, 0.95, 'reference', fontsize=10, backgroundcolor='#f5e3b3')
fig.text(0.47, 0.95, 'HeptA4T', fontsize=10, backgroundcolor='#f5e3b3')
fig.text(0.73, 0.95, 'SpaceG11T', fontsize=10, backgroundcolor='#f5e3b3')
fig.text(0.08, 0.95, '(A)', fontsize=10)
fig.text(0.38, 0.95, '(B)', fontsize=10)
fig.text(0.65, 0.95, '(C)', fontsize=10)
fig.text(0.08, 0.47, '(D)', fontsize=10)
fig.text(0.38, 0.47, '(E)', fontsize=10)
fig.text(0.65, 0.47, '(F)', fontsize=10)

# Define the mutant axes:
mut_ax = {'WT12rss':ppc_ax[0:2], '12HeptA4T':ppc_ax[2: 4], '12SpacG11T':ppc_ax[4:]}
compare_ax = {'WT12rss':ax0, '12HeptA4T':ax1, '12SpacG11T': ax2}
salt_ax = {'Ca': 0, 'Mg': 1}
colors = {'Ca': 'green', 'Mg':'purple'}

# Iterate through the mutants and salts and plot the data. 
for m, a in mut_ax.items():
    for s, ad in salt_ax.items():
        d = data[(data['mutant']==m) & (data['salt']==s)]
        x = np.sort((d['dwell_time_min'].values - (21/60)) * 60)
        y = np.arange(0, len(x)) / len(x)
        a[ad].step(x, y, 'k-', lw=0.75, label='data')
        a[ad].step([], [], '-', color=colors[s], label=r'$e^{-t/\tau} / \tau$', lw=0.75)
        a[ad].legend(loc='lower right', handlelength=0.75, fontsize=7)

        # Plot the PPC
        _ppc = ppc[(ppc['mutant']==m) & (ppc['salt']==s)]
        for g, d in _ppc.groupby(['sample_idx']):
            if g%1 == 0:
                a[ad].step(d['draws'] * 60, d['ecdf'], '-', lw=0.01, alpha=0.5,
            color=colors[s])
        
        # Plot the compared empirical CDFS
        if s == 'Ca':
            label = 'Ca$^2+$'
        else:
            label = 'Mg$^2+$'
        compare_ax[m].step(x, y,color=colors[s], label=label)

ax0.legend(loc='lower right', fontsize=8)
plt.subplots_adjust(hspace=0.45)
plt.savefig('CaMg_posterior_predictive.pdf', bbox_inches='tight', facecolor='white')
#%%
