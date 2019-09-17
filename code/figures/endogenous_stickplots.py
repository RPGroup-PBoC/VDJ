#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

# Load and trim the various datasets. 
cuts = pd.read_csv('../../data/pooled_cutting_probability.csv')
cuts = cuts[(cuts['salt']=='Mg') & (cuts['hmgb1']==80)]
dwell = pd.read_csv('../../data/compiled_dwell_times.csv')
dwell = dwell[(dwell['salt']=='Mg') & (dwell['hmgb1']==80)]
loops = pd.read_csv('../../data/compiled_loop_freq_bs.csv')
counts = loops[(loops['salt']=='Mg') & (loops['hmgb1']==80)]

#%% Compute the quartiles of dwell time and loops per bead
median_dwell = dwell.groupby('mutant')['dwell_time_min'].median().reset_index()
dwell_25 = dwell.groupby('mutant')['dwell_time_min'].quantile(0.25).reset_index()
dwell_75 = dwell.groupby('mutant')['dwell_time_min'].quantile(0.75).reset_index()
dwell_quartiles = pd.DataFrame()
dwell_quartiles['mutant'] = median_dwell['mutant']
dwell_quartiles['median_dwell_min'] = median_dwell['dwell_time_min']
dwell_quartiles['quartile_low'] = dwell_25['dwell_time_min']
dwell_quartiles['quartile_high'] = dwell_75['dwell_time_min']

cuts.rename(columns={'mode':'mean_p_cut'}, inplace=True)
dfs = [counts, dwell, cuts]
valid_dfs = [[], [], []]
for i, d in enumerate(dfs):
    d = d.copy()
    for g, _d in d.groupby(['mutant']):
        if ('Hept' not in g) & ('Spac' not in g) & ('Non' not in g):
            valid_dfs[i].append(_d) 

endo_counts = pd.concat(valid_dfs[0])
endo_dwell = pd.concat(valid_dfs[1])
endo_cuts = pd.concat(valid_dfs[2])
#%%
# Set up the figure canvas
fig, ax = plt.subplots(3, 1, figsize=(3.42, 5), sharex=True)

fig.text(-0.05, 0.88, '(A)', fontsize=8)
fig.text(-0.05, 0.61, '(B)', fontsize=8)
fig.text(-0.05, 0.36, '(C)', fontsize=8)
# Generate a mapping for position to mutant
muts = endo_counts['mutant'].unique()

map = {m:i for m, i in zip(muts, np.arange(1, 13))}
map = {'DFL1613':1, 'DFL161':2, 'V1-135':3, 'V9-120':4,
       'V10-96':5, 'V19-93':6, 'WT12rss':7, 'V4-55':8,
       'V5-43':9, 'V8-18':10, 'V6-17':11, 'V6-15':12}
muts = list(map.keys())
muts[6] = 'V4-57-1'
muts[0] = "DFL16.1-3'"
muts[1] = "DFL16.1-5'"

usage = {m:i for m, i in zip(muts, np.arange(1, 13))}
#usage = {3:' (6.5±0.6%)', 4:' (6.1±0.5%)',
#        5:' (6.5±0.5%)', 6:' (6.0±0.6%)',
#        7:' (0.22%)',8:' (2.7±0.6%)',
#        9:' (3.4±0.8%)', 10:' (0.0%)', 
#        11:' (6.8±0.5%)', 12:' (5.5%)'}
usage = {3:' (6.5%)', 4:' (6.1%)',
        5:' (6.5%)', 6:' (6.0%)',
        7:' (0.22%)',8:' (2.7%)',
        9:' (3.4%)', 10:' (0.0%)', 
        11:' (6.8%)', 12:' (5.5%)'}
use = list(usage.keys())
axlab = muts.copy()
for n in usage:
        axlab[n-1] = axlab[n-1] + usage[n]

# Format and add labels
for a in ax:
    a.yaxis.set_tick_params(labelsize=8)
ax[-1].set_xticks(np.arange(1, 13))
ax[-1].set_xticklabels(axlab, rotation=90, fontsize=8)
ax[0].set_ylabel('looping frequency', fontsize=8)
ax[1].set_ylabel('dwell time [min]', fontsize=8)
ax[2].set_ylabel('cutting probability', fontsize=8)

# set limits
for a in ax:
    a.set_xlim([0, len(muts) + 1])

ax[0].set_yticks([0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax[2].set_yticks([0, 0.2,  0.4, 0.6, 0.8, 1])
ax[1].set_yticks([0, 4, 8, 12, 16, 20])
ax[0].set_ylim([0, 0.65])
ax[1].set_ylim([0, 18])
ax[2].set_ylim([0, 1])

# Add stripes to separate the endogenous mutants
for a in ax:
    for i in np.arange(0, 13, 1):
        a.axvspan(i-0.25, i+0.25, color='white',
                 alpha=0.65, linewidth=0, zorder=-1)

# looping frequency
for g, d in endo_counts.groupby('mutant'):
    if g == 'WT12rss':
            face = 'dodgerblue' 
    else:
            face = 'w'
 
    ax[0].vlines(map[g], d['bs_low'], d['bs_high'], color='dodgerblue', lw=1)
    ax[0].plot(map[g], d['mean'], marker='o', markeredgecolor='dodgerblue', 
                markerfacecolor=face, ms=5)


"""# Boxplot
dwell_endog = pd.DataFrame()
for g, d in endo_dwell.groupby('mutant'):
        bp = ax[1].boxplot(d['dwell_time_min'],positions=[map[g]],
                        patch_artist=True)
        for patch in bp['boxes']:
                patch.set_facecolor('tomato')
                patch.set_alpha(0.5)"""

# Median dwell time
for g, d in endo_dwell.groupby('mutant'):
        if g == 'WT12rss':
                face = 'tomato' 
        else:
                face = 'w'

        quartiles = np.quantile(d['dwell_time_min'].values,
                                [0.5, 0.25, 0.75])

#        ax[1].vlines(map[g], 0, d['dwell_time_min'], color='tomato', lw=1)
        ax[1].plot(map[g], quartiles[0], marker='o', markeredgecolor='tomato', 
                markerfacecolor=face, ms=5)
        ax[1].vlines(map[g], quartiles[1], quartiles[2],
                        color='tomato', lw=2)

# Cutting probability
for g, d in endo_cuts.groupby('mutant'):
    if g == 'WT12rss':
            face = 'rebeccapurple' 
    else:
            face = 'w'
 
    ax[2].vlines(map[g], d['mean_p_cut']-d['std'], d['mean_p_cut']+d['std'], color='rebeccapurple', lw=1)
    ax[2].plot(map[g], d['mean_p_cut'], marker='o', markeredgecolor='rebeccapurple', 
                markerfacecolor=face, ms=5)


plt.savefig('./FigX_endogenous_properties.pdf', facecolor='white', bbox_inches='tight')
#%%
