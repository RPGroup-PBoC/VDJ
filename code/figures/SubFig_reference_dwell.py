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
fig,ax = plt.subplots(1, 1, figsize=(8.2,8))
bins=np.arange(0, 20, 1.0)
ax.hist(dwell['dwell_time_min'], color='#e28371', bins=bins)
ax.scatter(dwell['dwell_time_min'].median(), 38.5, color='slategrey',
            s=100)
ax.hlines(38.5, dwell['dwell_time_min'].quantile(0.25),
        dwell['dwell_time_min'].quantile(0.75), color='slategrey',
        lw=4, ls='-')

_ = ax.set_title('V4-57-1 (reference)', fontsize=24)
_ = ax.set_xlim([0, 20])
_ = ax.set_xticklabels([])
_ = ax.set_yticks(np.arange(0, 40, 10))
_ = ax.set_yticklabels([0, 10, 20, 30], fontsize=20)
_ = ax.set_ylabel('number of loops', fontsize=24)
_ = ax.set_xlabel('\ntime [min]', fontsize=24)
_ = ax.text(dwell['dwell_time_min'].quantile(0.25), 38, '25%',
            horizontalalignment='right', fontsize=14)
_ = ax.text(dwell['dwell_time_min'].quantile(0.75), 38, '75%',
            horizontalalignment='left', fontsize=14)
_ = ax.text(dwell['dwell_time_min'].quantile(0.5)+0.2, 39.2, '50%',
            horizontalalignment='center', fontsize=14)

fig.savefig('../../figures/SubFigB_reference_dwell_histogram.pdf', bbox_inches='tight',
            facecolor='white')

# %%
