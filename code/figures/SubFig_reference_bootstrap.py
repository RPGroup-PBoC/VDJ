"""
Representative Looping Frequency Bootstrap and 95% Confidence Interval
--------------------------------------------------------------------------------
Author: Soichi Hirokawa
Last Modified: January 7, 2020
License: MIT

Description
--------------------------------------------------------------------------------
This script generates the subfigure in the manuscript which shows a 
representative bootstrap replicate distribution and the 95% confidence interval.

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


# Upload V4-57-1 sequence looping dataset
data = pd.read_csv('../../data/compiled_looping_events.csv', comment='#')
data = data[(data['mutant']=='WT12rss') & (data['hmgb1']==80) & 
                (data['salt']=='Mg')]

percentiles = [2.5, 97.5]
col_names = ["bs_95_low", "bs_95_high"]
bs_reps = int(1E6)
bs_df = pd.DataFrame([])
sampling = np.random.choice(data['n_loops'].values,size=(len(data), bs_reps),
                            replace=True)
loop_freq = np.sum(sampling, axis=0) / len(data)
df_dict = {'mutant':'V4-57-1', 'salt':'Mg', 'hmgb1':80,
            'n_loops':data['n_loops'].sum(), 'n_beads':len(data),
            'loops_per_bead':data['n_loops'].sum() / len(data)}
computed_percentiles = np.percentile(loop_freq, percentiles)
for i,col in zip(computed_percentiles,col_names):
    df_dict[col] = i

bs_df = bs_df.append(df_dict, ignore_index=True)

# Form ECDFs
x = list(np.sort(loop_freq))
y = list(np.arange(0, bs_reps, 1) / bs_reps)
y_short = [-1, 2]
text_perc = '95%'

true_loops_val = y[x.index(bs_df['loops_per_bead'].values[0])]
#%%
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.set_ylim([-0.01, 1.01])
ax.step(x, y, color='dodgerblue', lw=2, zorder=10)
ax.axvline(bs_df['loops_per_bead'].values[0], 0, true_loops_val, 
            color='slategrey', ls='--', alpha=0.4, lw=2)
for n in range(2):
    ax.hlines(percentiles[n]/100, 0, bs_df[col_names[n]].values[0],linestyle='--',
                color='slategrey', alpha=0.4) 
    ax.vlines(bs_df[col_names[n]].values[0], percentiles[n]/100, 
                true_loops_val, linestyle='--', color='slategrey', 
                alpha=0.4)
ax.scatter(bs_df['loops_per_bead'], true_loops_val, color='slategrey',
            s=50, alpha=0.7)
ax.hlines(true_loops_val, bs_df[col_names[0]].values[0],
        bs_df[col_names[1]].values[0], alpha=0.7, 
        ls='-', color='slategrey', lw=3)

_ = ax.set_xlim([loop_freq.min(),loop_freq.max()])
_ = ax.set_xlabel('looping frequency', fontsize=12)
_ = ax.set_ylabel('ECDF', fontsize=12)
_ = ax.set_xticks([0.15, 0.20, 0.250, 0.30])
_ = ax.set_xticklabels([0.15, 0.20, 0.25, 0.30])

fig.savefig('../../figures/SubFigXB_reference_bootstrap.pdf',
            bbox_inches='tight', facecolor='white')
# %%
