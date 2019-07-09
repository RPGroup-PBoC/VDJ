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

cut_data = pd.read_csv('../../data/pooled_cutting_probability.csv')
cut_data = cut_data[(cut_data['hmgb1'] == 80) & (cut_data['salt']=='Mg')]


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

colors = {'A':'#E10C00', 'T':'#38C2F2', 'C':'#278C00', 'G':'#5919FF'}
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
#plt.savefig('./loop_frequency_stickplot.pdf', facecolor='white')
#%%
# Obtain information on endogenous sequences

# Rename some of the mutant names to correct endogenous names
endog_names = {'WT12rss' : 'V4-57-1 (ref)',
                '12SpacC1A' : 'V4-55',
                'DFL161' : 'DFL 16.1-5\'',
                'DFL1613' : 'DFL 16.1-3\''}

# Isolate only the point mutants for now
endog = cut_data[(~cut_data['mutant'].str.startswith('12')) | (cut_data['mutant']=='12SpacC1A')].copy()
endog_counts = counts[(~counts['mutant'].str.startswith('12')) | (counts['mutant']=='12SpacC1A')].copy()

endog = endog.replace({'mutant' : endog_names})
endog_counts = endog_counts.replace({'mutant' : endog_names})

endog_ordering_y = {'V6-15' : 0.0,
                    'V6-17' : 1.0,
                    'V8-18' : 2.0,
                    'V5-43' : 3.0,
                    'V4-55' : 4.0,
                    'V4-57-1 (ref)' : 5.0,
                    'V19-93' : 6.0,
                    'V10-96' : 7.0,
                    'V9-120' : 8.0,
                    'V1-135' : 9.0,
                    'DFL 16.1-5\'' : 10.0,
                    'DFL 16.1-3\'' : 11.0,
                    }

fig, ax = plt.subplots(1, 2, figsize=(6, 6))

for mut, mut_info in endog_counts.groupby('mutant'):
    if mut_info['rel_diff'].values[0] < 0:
        ax[0].plot(mut_info['rel_diff'], endog_ordering_y[mut], marker='o',
                color='#E10C00', ms=4.5, linestyle='None')
        ax[0].hlines(endog_ordering_y[mut], 0, mut_info['rel_diff'],  color='#E10C00',
                linewidth=2, label='__nolegend__')
    elif mut_info['rel_diff'].values[0] > 0:
        ax[0].plot(mut_info['rel_diff'], endog_ordering_y[mut], marker='o',
                color='#38C2F2', ms=4.5, linestyle='None')
        ax[0].hlines(endog_ordering_y[mut], 0, mut_info['rel_diff'],  color='#38C2F2',
                linewidth=2, label='__nolegend__')
    else:
        ax[0].plot(mut_info['rel_diff'], endog_ordering_y[mut], marker='o',
                color='k', ms=4.5, linestyle='None')
        ax[0].hlines(endog_ordering_y[mut], 0, mut_info['rel_diff'],  color='k',
                linewidth=2, label='__nolegend__')

ref = endog[endog['mutant']=='V4-57-1 (ref)']
ref_mode = ref['mode'].values[0]

for mut, mut_info in endog.groupby('mutant'):
    if mut_info['mode'].values[0] - ref_mode < 0:
        ax[1].plot(mut_info['mode'] - ref_mode, endog_ordering_y[mut], marker='o',
                color='#E10C00', ms=4.5, linestyle='None')
        ax[1].hlines(endog_ordering_y[mut], 0, mut_info['mode'] - ref_mode,  color='#E10C00',
                linewidth=2, label='__nolegend__')
    elif mut_info['mode'].values[0] - ref_mode > 0:
        ax[1].plot(mut_info['mode'] - ref_mode, endog_ordering_y[mut], marker='o',
                color='#38C2F2', ms=4.5, linestyle='None')
        ax[1].hlines(endog_ordering_y[mut], 0, mut_info['mode'] - ref_mode,  color='#38C2F2',
                linewidth=2, label='__nolegend__')
    else:
        ax[1].plot(mut_info['mode'] - ref_mode, endog_ordering_y[mut], marker='o',
                color='k', ms=4.5, linestyle='None')
        ax[1].hlines(endog_ordering_y[mut], 0, mut_info['mode'] - ref_mode,  color='k',
                linewidth=2, label='__nolegend__')

for a in ax:
        _ = a.set_yticks(np.arange(0, 12))
        a.set_ylim([-0.5, 11.5])
        a.vlines(0, 0, len(endog_ordering_y), color='k', linestyle=':', zorder=-1)
        for i in range(0, 12, 2):
                a.axhspan(i-0.5, i+0.5, color='w', linewidth=0, zorder=-2)

_ = ax[0].set_yticklabels(list(endog_ordering_y), fontsize=12)
_ = ax[1].set_yticklabels([])
ax[0].set_xlim([-0.25, 0.25])
ax[1].set_xlim([-0.45, 0.45])

ax[0].set_xlabel('change in\nloop frequency', fontsize=14)
ax[1].set_xlabel('change in cut\nprobability', fontsize=14)

plt.savefig('./endog_stickplot.pdf', facecolor='white', bbox_inches='tight')

#%%
