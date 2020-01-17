"""
Dwell Time Distribution and Quartiles Subfigure
--------------------------------------------------------------------------------
Author: Soichi Hirokawa
Last Modified: January 7, 2020
License: MIT

Description
--------------------------------------------------------------------------------
This script generates the subfigure for the dwell time distribution histogram
and how the quartiles are identified for subfigure panel C

Notes
--------------------------------------------------------------------------------
This script is designed to be executed from the `code/figures` directory and uses 
a relative path to load the necessary CSV files.
"""
import numpy as np
import pandas as pd
import vdj.io
import vdj.viz
import matplotlib.pyplot as plt
import matplotlib.patches as patches
vdj.viz.plotting_style()

# Load the dwell times
dwell = pd.read_csv('../../data/compiled_dwell_times.csv', comment='#')
dwell = dwell[(dwell['salt']=='Mg') & (dwell['hmgb1']==80) & (dwell['mutant']=='WT12rss')]

#%%
fig,ax = plt.subplots(1, 1, figsize=(9,3))
bins=np.arange(0, 20, 1.0)
ax.hist(dwell['dwell_time_min'], color='#e28371', bins=bins)
ax.plot(dwell['dwell_time_min'].median(), 15, color='dodgerblue', lw=1.25,
            ms=25, zorder=10, marker='o', markerfacecolor='white')
ax.hlines(15, dwell['dwell_time_min'].quantile(0.25),
        dwell['dwell_time_min'].quantile(0.75), color='dodgerblue',
        lw=5, ls='-', zorder=10)
ax.text(dwell['dwell_time_min'].median(), 14.6, 'N', fontsize=18, zorder=20,
        color='dodgerblue', horizontalalignment='center',
        verticalalignment='center')

_ = ax.set_xticks(np.arange(0, 30, 5))
_ = ax.set_xticklabels(np.arange(0, 30, 5), fontsize=18)
_ = ax.set_xlim([0, 20])
_ = ax.set_yticks(np.arange(0, 40, 10))
_ = ax.set_yticklabels(np.arange(0, 40, 10), fontsize=18)
_ = ax.set_ylabel('counts', fontsize=28)
_ = ax.set_xlabel('time [min]', fontsize=28)

fig.savefig('../../figures/SubFigB_reference_dwell_histogram.pdf', bbox_inches='tight',
            facecolor='white')

# %%
