# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

# Load the dwell time measurements and restrict 
data = pd.read_csv('../../data/compiled_dwell_times.csv', comment='#')
data = data[(data['mutant']=='WT12rss') | (data['mutant']=='12SpacG11T') | 
            (data['mutant']=='12HeptA4T')].copy()

# Load the sampling statistics and restrict
stats = pd.read_csv('../../data/expon_waiting_time_posterior_summary.csv', comment='#')
stats = stats[(stats['mutant']=='WT12rss') | (stats['mutant']=='12SpacG11T') | 
              (stats['mutant']=='12HeptA4T')]

#%%
# Set up the figure canvas using gridspec
fig = plt.figure(figsize=(7, 4))
gs = GridSpec(6, 6)
ax1 = fig.add_subplot(gs[0:2, 0])
ax2 = fig.add_subplot(gs[0:2, 1])
ax3 = fig.add_subplot(gs[0:2, 2])
ax4 = fig.add_subplot(gs[0:2, 3])
ax5 = fig.add_subplot(gs[0:2, 4])
ax6 = fig.add_subplot(gs[0:2, 5])
ax7 = fig.add_subplot(gs[3:,:2])
ax8 = fig.add_subplot(gs[3:,2:4])
ax9 = fig.add_subplot(gs[3:,4:])

# Assemble the axes for easy iteration
ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]

# Define the axes by mutant
ax_map = {'WT12rss':[ax1, ax2, ax7], '12HeptA4T':[ax3, ax4, ax8],
          '12SpacG11T': [ax5, ax6, ax9]}

# Format the labels 
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('dwell time\n[min]', fontsize=8)
    a.set_xscale('log')
    a.set_ylim([0, 1])
    a.set_xlim([0.35, 80])
    if (i != 0) & (i != 6):
        a.set_yticks([])
for i, a in enumerate(ax[:6]):
    if i%2 != 0:
        title = 'Mg$^{2+}$'
    else:
        title = 'Ca$^{2+}$'
    a.set_title(title, fontsize=8)
ax[0].set_ylabel('cumulative\ndistribution', fontsize=8)
ax[6].set_ylabel('cumulative distribution', fontsize=8)

# Add panel labels and titles.
fig.text(0.2, 0.95, 'V4-57-1 (ref)', fontsize=9)
fig.text(0.45, 0.95, 'Heptamer A4T', fontsize=9)
fig.text(0.72, 0.95, 'Spacer G11T', fontsize=9)
fig.text(0.05, 0.95, '(A)', fontsize=9)
fig.text(0.4, 0.95, '(B)', fontsize=9)
fig.text(0.65, 0.95, '(C)', fontsize=9)

# Plot the empirical CDFs
DEADFILTER = 21/60
time_range = np.linspace(0, 80, 500)  - DEADFILTER

for g, d in data.groupby(['mutant']):
    # Parse the axis mapping
    a = ax_map[g]

    # Separate by salt
    ca = d[d['salt']=='Ca']
    mg = d[d['salt']=='Mg']

    # Compute the distributions    
    ca_x = np.sort(ca['dwell_time_min'].values) - DEADFILTER
    ca_y = np.arange(0, len(ca), 1) / len(ca)
    mg_x = np.sort(mg['dwell_time_min'].values) - DEADFILTER
    mg_y = np.arange(0, len(mg), 1) / len(mg)

    # Parse the statistics 
    ca_stats = stats[(stats['mutant']==g) & (stats['salt']=='Ca')]
    mg_stats = stats[(stats['mutant']==g) & (stats['salt']=='Mg')]

    # Compute the exponential dists 
    ca_low = 1 - np.exp(-(time_range - DEADFILTER)/ca_stats['hpd_min'].values[0])
    ca_high = 1 - np.exp(-(time_range - DEADFILTER)/ca_stats['hpd_max'].values[0])
    mg_low = 1 - np.exp(-(time_range - DEADFILTER)/mg_stats['hpd_min'].values[0])
    mg_high = 1 - np.exp(-(time_range - DEADFILTER)/mg_stats['hpd_max'].values[0])

    # Plot the fits. 
    a[0].fill_between(time_range, ca_low, ca_high, color='green', alpha=0.45,
                    label='fit')
    a[1].fill_between(time_range, mg_low, mg_high, color='rebeccapurple', alpha=0.45,
                    label='fit')

    # Plot the cdfs
    a[0].step(ca_x, ca_y, 'k-', lw=1, label='data')
    a[1].step(mg_x, mg_y, 'k-', lw=1, label='data')
    a[2].step(ca_x, ca_y, color='green', lw=1, label='Ca$^{2+}$ data')
    a[2].step(mg_x, mg_y, color='rebeccapurple', lw=1, label='Mg$^{2+}$ data')

ax[6].legend(fontsize=8, handlelength=0.75)
ax[0].legend(fontsize=7, handlelength=0.75)
plt.savefig('../../figures/FigX_CaMg_distributions.pdf', bbox_inches='tight', face_color='white')

#%%


#%%
