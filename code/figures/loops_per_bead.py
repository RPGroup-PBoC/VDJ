#-*-coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as lines
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
fig, ax = plt.subplots(1, 1, figsize=(7, 4))

colors = {'A':'#E10C00', 'T':'#38C2F2', 'C':'#36C200', 'G':'#5919FF'}
shift = {'A':-0.26, 'T':0.26, 'C':-0.13, 'G':0.13}

for g, d in points.groupby(['base']):
    ax.plot(d['pos'] + shift[g] + 1, d['rel_diff'], marker='o', color=colors[g], lw=1, 
                ms=4.5, linestyle='none', label=g)
    ax.vlines(d['pos'] + shift[g] + 1, 0, d['rel_diff'], color=colors[g], lw=2, label='__nolegend__')

_ = ax.set_xticks(np.arange(1, 29))
_ = ax.set_xticklabels(list(ref_seq))
ax.set_xlim([0.5, 28.5])

line1 = lines.Line2D([7.5, 7.5], [-0.35, 0.-0.31], clip_on=False, alpha=1,
                    linewidth=1, color='k')
line2 = lines.Line2D([19.5, 19.5], [-0.35, 0.-0.31], clip_on=False, alpha=1,
                    linewidth=1, color='k')
ax.add_line(line1)
ax.add_line(line2)


ax.legend(fontsize=8, ncol=5)
ax.set_xlabel('reference sequence', fontsize=12)
ax.hlines(0, 0, 29, color='k', linestyle=':')
ax.set_ylim([-0.3, 0.3])
ax.set_ylabel('change in\nloop frequency', fontsize=12)
for i in range(1, 29, 2):
    ax.axvspan(i-0.5, i+0.5, color='w', linewidth=0, zorder=-1)
ax.vlines(0.5, -0.4, 0.4, color='#f5e3b3', linewidth=4, zorder=0)
plt.savefig('./loop_frequency_stickplot.pdf', facecolor='white')
#%%





#%%
