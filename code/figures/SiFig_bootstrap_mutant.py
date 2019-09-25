"""
Representative Bootstrap Replicates of Looping Frequency
--------------------------------------------------------------------------------
Author: Soichi Hirokawa
Last Modified: September 25, 2019
License: MIT

Description
--------------------------------------------------------------------------------
This script generates SI Figure 2 which shows a representative bootstrap
replicate distribution and confidence intervals.

Notes
--------------------------------------------------------------------------------
This script is designed to be executed from the `code/figures` directory and uses a relative path to load the necessary CSV files.
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

percentiles = [(2.5, 97.5), (87.5, 12.5), (25, 75), (37.5, 62.5), (45, 55), 
               (47.5, 52.5)]
col_names = [("bs_95_low", "bs_95_high"), ("bs_75_low", "bs_75_high"),
             ("bs_50_low", "bs_50_high"),  ("bs_25_low", "bs_25_high"),
             ("bs_10_low", "bs_10_high"), ("bs_5_low", "bs_5_high")]
bs_reps = int(1E6)
bs_df = pd.DataFrame([])
sampling = np.random.choice(data['n_loops'].values,size=(len(data), bs_reps),
                            replace=True)
loop_freq = np.sum(sampling, axis=0) / len(data)
df_dict = {'mutant':'V4-57-1', 'salt':'Mg', 'hmgb1':80,
            'n_loops':data['n_loops'].sum(), 'n_beads':len(data),
            'loops_per_bead':data['n_loops'].sum() / len(data)}
for perc, col in zip(percentiles, col_names):
    computed_percentiles = np.percentile(loop_freq, perc)
    for i in range(2):
        df_dict[col[i]] = computed_percentiles[i]

bs_df = bs_df.append(df_dict, ignore_index=True)

# Form ECDFs
x = list(np.sort(loop_freq))
y = list(np.arange(0, bs_reps, 1) / bs_reps)
y_short = [-1, 2]
heights = [0.05, 0.1, 0.15, 0.20, 0.25, 0.30]
text_perc = ['5%', '10%', '25%', '50%', '75%', '95%']

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.set_ylim([-0.01, 1.01])
ax.plot(x, y, color='dodgerblue', lw=2)
ax.axvline(bs_df['loops_per_bead'].values[0], 0, 1.0, color='grey', lw=2)
for percentile in col_names:
    ax.fill_betweenx(y_short, bs_df[percentile[0]].values[0], bs_df[[percentile[1]]].values[0],
                    color='grey', alpha=0.2)
ax.set_xlabel('looping frequency', fontsize=12)
ax.set_ylabel('ECDF', fontsize=12)
for n in range(len(heights)):
    ax.add_patch(patches.Rectangle((0.295,0.05),0.025,heights[n],
                facecolor='grey', alpha=0.2))
    ax.text(0.32, 0.085 + (n - 0.5) * heights[0], text_perc[n])
plt.savefig('../../figures/SiFig_bootstrap_reference.pdf', facecolor='white',
            bbox_inches='tight')

