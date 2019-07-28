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
loops = pd.read_csv('../../data/compiled_looping_events.csv')
loops = loops[(loops['salt']=='Mg') & (loops['hmgb1']==80)]

#%% Compute the median dwell time and loops per bead
median_dwell = dwell.groupby('mutant')['dwell_time_min'].median().reset_index()
counts = loops.groupby(['mutant'])[['n_loops']].agg(('sum', 'count')).reset_index()
counts['loops_per_bead'] = counts['n_loops']['sum'] / counts['n_loops']['count']
cuts.rename(columns={'mode':'mean_p_cut'}, inplace=True)
dfs = [counts, median_dwell, cuts]
valid_dfs = [[], [], []]
for i, d in enumerate(dfs):
    d = d.copy()
    d[d['mutant']=='12SpacC1A', 'mutant'] = 'V4-55'
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

# Format and add labels
for a in ax:
    a.yaxis.set_tick_params(labelsize=8)
ax[-1].set_xticks(np.arange(1, 13))
ax[-1].set_xticklabels(muts, rotation=90, fontsize=8)
ax[0].set_ylabel('looping frequency', fontsize=8)
ax[1].set_ylabel('median\ndwell time [min]', fontsize=8)
ax[2].set_ylabel('cutting probability', fontsize=8)

# set limits
for a in ax:
    a.set_xlim([0, len(muts) + 1])

ax[0].set_yticks([0, 0.2,  0.4, 0.6])
ax[2].set_yticks([0, 0.2,  0.4, 0.6, 0.8, 1])
ax[1].set_yticks([0, 1, 2, 3, 4, 5])
ax[0].set_ylim([0, 0.75])
ax[1].set_ylim([0, 5])
ax[2].set_ylim([0, 1])

# Add stripes to separate the endogenous mutants
for a in ax:
    for i in np.arange(0, 13, 1):
        a.axvspan(i-0.25, i+0.25, color='white',
                 alpha=0.65, linewidth=0, zorder=-1)

# looping frequency
for g, d in endo_counts.groupby('mutant'):
    ax[0].vlines(map[g], 0, d['loops_per_bead'], color='dodgerblue', lw=1)
    ax[0].plot(map[g], d['loops_per_bead'], marker='o', markeredgecolor='dodgerblue', 
                markerfacecolor='white', ms=5)

# Median dwell time
for g, d in endo_dwell.groupby('mutant'):
    ax[1].vlines(map[g], 0, d['dwell_time_min'], color='tomato', lw=1)
    ax[1].plot(map[g], d['dwell_time_min'], marker='o', markeredgecolor='tomato', 
                markerfacecolor='white', ms=5)

# Cutting probability
for g, d in endo_cuts.groupby('mutant'):
    ax[2].vlines(map[g], 0, d['mean_p_cut'], color='rebeccapurple', lw=1)
    ax[2].plot(map[g], d['mean_p_cut'], marker='o', markeredgecolor='rebeccapurple', 
                markerfacecolor='white', ms=5)

plt.savefig('./FigX_endogenous_properties.pdf', facecolor='white', bbox_inches='tight')
#%%
