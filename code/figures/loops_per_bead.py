#-*-coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.io 
import vdj.viz
vdj.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/compiled_looping_events.csv')
data = data[(data['salt']=='Mg') & (data['hmgb1']==80)]


#%%
counts = data.groupby(['mutant'])[['n_loops']].agg(('sum', 'count')).reset_index()
counts['loops_per_bead'] = counts['n_loops']['sum'] / counts['n_loops']['count']

# Get the reference seq
ref = vdj.io.endogenous_seqs()['WT12rss']
ref_seq = ref[0]
ref_idx = ref[1]

# Include the mutant id information
wt_val = counts[counts['mutant']=='WT12rss']['loops_per_bead'].values[0]
for m in counts['mutant'].unique():
    seq = vdj.io.mutation_parser(m)
    counts.loc[counts['mutant']==m, 'n_muts'] = seq['n_muts']
    
    # Find the x and mutation identity
    loc = np.argmax(ref_idx != seq['seq_idx'])
    mut = seq['seq'][loc]
    counts.loc[counts['mutant']==m, 'pos'] = loc
    counts.loc[counts['mutant']==m, 'base'] = mut


for m in counts['mutant'].unique():
    _d = counts[counts['mutant']==m]
    if _d['loops_per_bead'].values[0] < wt_val:
        val = 1 - _d['loops_per_bead'] / wt_val
    else:
        val = _d['loops_per_bead'] / wt_val - 1
    counts.loc[counts['mutant']==m, 'rel_diff'] = _d['loops_per_bead'].values[0] - wt_val 

# Keep the single point mutants
points = counts[counts['n_muts'] == 1].copy()

#%%
bar_width = 0.75
fig, ax = plt.subplots(2, 1, figsize=(7, 4))

colors = {'A':'tomato', 'T':'dodgerblue', 'C':'grey', 'G':'purple'}
shift = {'A':-0.25, 'T':0.15, 'C':-0.15, 'G':0.25}

for g, d in points.groupby(['base']):
    ax[1].plot(d['pos'] + shift[g] + 1, d['rel_diff'], marker='o', color=colors[g], lw=1, 
                ms=3, linestyle='none', label=g, zorder=order)
    ax[1].vlines(d['pos'] + shift[g] + 1, 0, d['rel_diff'], color=colors[g], lw=1, label='__nolegend__')

for g, d in points.groupby(['base']):
    ax[0].plot(d['pos'] + shift[g] + 1, d['loops_per_bead'], marker='o', color=colors[g], lw=1, 
                ms=3, linestyle='none', label=g, zorder=order)
    ax[0].vlines(d['pos'] + shift[g] + 1, 0, d['loops_per_bead'], color=colors[g], lw=1, label='__nolegend__')

ax[0].hlines(wt_val, 0, 29, linestyle=':', color='k', label='WT')
for a in ax:
    _ = a.set_xticks(np.arange(1, 29))
    _ = a.set_xticklabels(list(ref_seq))
    a.set_xlim([0, 29])
    a.legend(fontsize=8, ncol=5)
    a.set_xlabel('reference sequence')
ax[1].hlines(0, 0, 29, color='k', linestyle=':')
ax[1].set_ylim([-0.3, 0.3])
ax[0].set_ylim([-0.01, 0.6])
ax[1].set_ylabel('change in\nloop frequency')
ax[0].set_ylabel('loop frequency')
for i in range(1, 29, 2):
    ax[0].vlines(i, -0.4, 0.8, color='w', linewidth=12, zorder=-1)
    ax[1].vlines(i, -0.4, 0.8, color='w', linewidth=12, zorder=-1)

plt.savefig('./loop_frequency_stickplot.pdf')
#%%





#%%
