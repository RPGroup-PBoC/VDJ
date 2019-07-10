#-*-coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import seaborn as sns
import vdj.io 
import vdj.viz
vdj.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/compiled_looping_events.csv')
data = data[(data['salt']=='Mg') & (data['hmgb1']==80)]

cut_data = pd.read_csv('../../data/pooled_cutting_probability.csv')
cut_data = cut_data[(cut_data['hmgb1'] == 80) & (cut_data['salt']=='Mg')]

cut_posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv')
cut_posts = cut_posts[(cut_posts['hmgb1']==80) & (cut_posts['salt']=='Mg')]

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
#%%
counts = data.groupby(['mutant'])[['n_loops']].agg(('sum', 'count')).reset_index()
counts['loops_per_bead'] = counts['n_loops']['sum'] / counts['n_loops']['count']

# Get the reference seq
ref = vdj.io.endogenous_seqs()['WT12rss']
ref_seq = ref[0]
ref_idx = ref[1]

# Include the mutant id information
wt_val = counts[counts['mutant']=='WT12rss']['loops_per_bead'].values[0]
wt_cut = cut_data[cut_data['mutant']=='WT12rss']['mode'].values[0]
wt_std = cut_data[cut_data['mutant']=='WT12rss']['std'].values[0]
for m in counts['mutant'].unique():
    seq = vdj.io.mutation_parser(m)
    counts.loc[counts['mutant']==m, 'n_muts'] = seq['n_muts']
    cut_data.loc[cut_data['mutant']==m, 'n_muts'] = seq['n_muts']

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
points_cut = cut_data[cut_data['n_muts'] == 1].copy()

for m in points_cut['mutant'].unique():
        seq = vdj.io.mutation_parser(m)
        loc = np.argmax(ref_idx != seq['seq_idx'])
        mut = seq['seq'][loc]
        points_cut.loc[points_cut['mutant']==m, 'pos'] = loc
        points_cut.loc[points_cut['mutant']==m, 'base'] = mut
        _d = points_cut[points_cut['mutant']==m]
        points_cut.loc[points_cut['mutant']==m, 'rel_diff'] = _d['mode'].values[0] - wt_cut

posterior_list = ['12HeptC3G', '12HeptC3T', '12SpacC4G', '12NonA4C', 'WT12rss', '12SpacG10T']
post_colors = {'12HeptC3G' : '#E10C00', 
                '12HeptC3T' : '#BF3030', 
                '12SpacC4G' : '#A24F59', 
                '12NonA4C' : '#7D778E', 
                'WT12rss' : '#599DC1', 
                '12SpacG10T' : '#38C2F2'}
post_hatch = {'12HeptC3G' : None, 
                '12HeptC3T' : None, 
                '12SpacC4G' : None, 
                '12NonA4C' : None, 
                'WT12rss' : None, 
                '12SpacG10T' : None}
plot_offset = dict(zip(posterior_list[::-1], np.arange(0.0, 0.3, 0.3/len(posterior_list))))
#%%
bar_width = 0.75
fig, ax = plt.subplots(3, 1, figsize=(8, 7))
plt.subplots_adjust(hspace=0.4)

# Isolate cut information with few beads for points_cut
low_cut = points_cut[points_cut['n_beads']<=20].copy()
high_cut = points_cut[points_cut['n_beads']>=20].copy()

colors = {'A':'#E10C00', 'T':'#38C2F2', 'C':'#278C00', 'G':'#5919FF'}
shift = {'A':-0.26, 'T':0.26, 'C':-0.13, 'G':0.13}

for g, d in points.groupby(['base']):
    ax[0].plot(d['pos'] + shift[g] + 1, d['rel_diff'], marker='o', color=colors[g], lw=1, 
                ms=4.5, linestyle='none', label=g)
    ax[0].vlines(d['pos'] + shift[g] + 1, 0, d['rel_diff'], color=colors[g], lw=2, label='__nolegend__')

for g, d in high_cut.groupby(['base']):
        ax[1].plot(d['pos'] + shift[g] + 1, d['mode']-wt_cut, marker='o', color=colors[g], lw=1,
                ms=4.5, linestyle='none', label=g)
        ax[1].vlines(d['pos'] + shift[g] + 1, 0, d['mode']-wt_cut, color=colors[g], lw=1,
                label='__nolegend__', linewidth=2)

# Plot HeptT6A:
for g, d in low_cut.groupby(['base']):
        ax[1].plot(d['pos'] + shift[g] + 1, d['mode']-wt_cut, marker='o', color=colors[g], lw=1,
                        ms=4.5, linestyle='none', label=g, alpha=0.3)
        ax[1].vlines(d['pos'] + shift[g] + 1, 0, d['mode']-wt_cut, color=colors[g], lw=1,
                        label='__nolegend__', alpha=0.3, linewidth=2)

line1 = lines.Line2D([7.5, 7.5], [-0.84, -0.72], clip_on=False, alpha=1,
                    linewidth=1, color='k')
line2 = lines.Line2D([19.5, 19.5], [-0.84, -0.72], clip_on=False, alpha=1,
                    linewidth=1, color='k')

for n in range(0,2):
        _ = ax[n].set_xticks(np.arange(1, 29))
        ax[n].set_xlim([0.5, 28.5])
        ax[n].vlines(0.5, -0.65, 1.0, color='#f5e3b3', linewidth=4, zorder=0)
        ax[n].hlines(0, 0, 29, color='k', linestyle=':')
        for i in range(1, 29, 2):
                ax[n].axvspan(i-0.5, i+0.5, color='w', linewidth=0, zorder=-1)

_ = ax[0].set_xticklabels([])
_ = ax[1].set_xticklabels(list(ref_seq))
ax[1].add_line(line1)
ax[1].add_line(line2)

ax[0].legend(fontsize=8, ncol=5)
ax[1].set_xlabel('reference sequence', fontsize=12)
ax[0].set_xlabel(None)
ax[0].set_ylim([-0.3, 0.3])
ax[1].set_ylim([-0.65, 0.65])
ax[0].set_ylabel('change in\nloop frequency', fontsize=12)
ax[1].set_ylabel('change in cut\nprobability', fontsize=12)

df_post = cut_posts.loc[cut_posts['mutant'].isin(posterior_list)]

sort_index = dict(zip(posterior_list, range(len(posterior_list))))
df_post['rank_index'] = df_post['mutant'].map(sort_index)
df_post.sort_values(['rank_index', 'probability'], ascending=True, inplace=True)
df_post.drop('rank_index', 1, inplace=True)

for mut, mut_posts in df_post.groupby('mutant'):
        ax[2].fill_between(mut_posts['probability'], plot_offset[mut],
                        mut_posts['posterior'] + plot_offset[mut],
                        color=post_colors[mut], alpha=1.0)
        ax[2].plot(mut_posts['probability'], mut_posts['posterior'] + plot_offset[mut],
                        color='white')
        ax[2].axhline(plot_offset[mut], 0, 1.0, color=post_colors[mut], alpha=1.0)
        if mut=='WT12rss':
                ax[2].text(1.00, 0.02 + plot_offset[mut], 'reference', 
                        fontsize=12, color='k', ha="right", va="center")
        else:
                ax[2].text(1.00, 0.02 + plot_offset[mut], mut, 
                        fontsize=12, color='k', ha="right", va="center")

ax[2].set_xlabel('probability of cut')
ax[2].set_ylim([-0.01, 0.4])
ax[2].set_xlim([0.0, 1.0])
ax[2].set_yticklabels([])
ax[2].set_ylabel('posterior distribution')

plt.savefig('./point_mutation_stickplot.pdf', facecolor='white', bbox_inches='tight')

#%%
fig, ax = plt.subplots(1, 1, figsize=(7,4))
for g, d in high_cut.groupby(['base']):
        ax.errorbar(d['pos'] + shift[g] + 1, d['mode'], yerr=d['std'], marker='o', color=colors[g], lw=1,
                ms=4.5, linestyle='none', label=g)

# Plot HeptT6A:
for g, d in low_cut.groupby(['base']):
        ax.errorbar(d['pos'] + shift[g] + 1, d['mode'], yerr=d['std'], marker='o', color=colors[g], lw=1,
                        ms=4.5, linestyle='none', alpha=0.3)

line1 = lines.Line2D([7.5, 7.5], [-0.08, -0.02], clip_on=False, alpha=1,
                    linewidth=1, color='k')
line2 = lines.Line2D([19.5, 19.5], [-0.08, -0.02], clip_on=False, alpha=1,
                    linewidth=1, color='k')

_ = ax.set_xticks(np.arange(1, 29))
ax.set_xlim([0.5, 28.5])
ax.vlines(0.5, -0.65, 1.0, color='#f5e3b3', linewidth=4, zorder=0)
for i in range(1, 29, 2):
        ax.axvspan(i-0.5, i+0.5, color='w', linewidth=0, zorder=-1)

_ = ax.set_xticklabels(list(ref_seq))
ax.add_line(line1)
ax.add_line(line2)
ax.axhspan(wt_cut-wt_std, wt_cut+wt_std, 0, 29, color='gray', 
        linestyle=':', alpha=0.4, label=r'$p_{cut}^{ref} \pm \sigma$')

ax.legend(fontsize=8, ncol=5)
ax.set_xlabel('reference sequence', fontsize=12)
ax.set_ylim([0.0, 1.0])
ax.set_ylabel('change in cut\nprobability', fontsize=12)

plt.savefig('./p_cut_point_SI.pdf', facecolor='white')
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

fig, ax = plt.subplots(1, 2, figsize=(6, 6))
plt.subplots_adjust(wspace=0.1)

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
ref_std = ref['std'].values[0]

for mut, mut_info in endog.groupby('mutant'):
        if mut=='V8-18':
                continue
        if mut_info['mode'].values[0] - ref_mode < 0:
                ax[1].plot(mut_info['mode'], endog_ordering_y[mut], marker='o',
                        color='#E10C00', ms=4.5, linestyle='None')
                ax[1].errorbar(mut_info['mode'], endog_ordering_y[mut], xerr=mut_info['std'], color='#E10C00',
                        linewidth=2, label='__nolegend__')
        elif mut_info['mode'].values[0] - ref_mode > 0:
                ax[1].plot(mut_info['mode'], endog_ordering_y[mut], marker='o',
                        color='#38C2F2', ms=4.5, linestyle='None')
                ax[1].errorbar(mut_info['mode'], endog_ordering_y[mut], xerr=mut_info['std'], color='#38C2F2',
                        linewidth=2, label='__nolegend__')
        else:
                ax[1].plot(mut_info['mode'], endog_ordering_y[mut], marker='o',
                        color='k', ms=4.5, linestyle='None')
                ax[1].errorbar(mut_info['mode'], endog_ordering_y[mut], xerr=mut_info['std'], color='k',
                        linewidth=2, label='__nolegend__')

ax[0].hlines(-0.5, -0.3, 0.3, color='#f5e3b3', linewidth=4, zorder=0)
ax[1].hlines(-0.5, -0.45, 0.45, color='#f5e3b3', linewidth=4, zorder=0)

for a in ax:
        _ = a.set_yticks(np.arange(0, 12))
        a.set_ylim([-0.5, 11.5])
        for i in range(0, 12, 2):
                a.axhspan(i-0.5, i+0.5, color='w', linewidth=0, zorder=-2)
ax[0].vlines(0, 0, len(endog_ordering_y), color='k', linestyle=':', zorder=-1)
ax[1].axvspan(ref_mode-ref_std, ref_mode+ref_std, 0, len(endog_ordering_y), 
        color='gray', linestyle=':', zorder=-1, alpha=0.4)
_ = ax[0].set_yticklabels(list(endog_ordering_y), fontsize=12)
_ = ax[1].set_yticklabels([])
ax[0].set_xlim([-0.25, 0.25])
ax[1].set_xlim([0.0, 1.0])

ax[0].set_xlabel('change in\nloop frequency', fontsize=14)
ax[1].set_xlabel('change in cut\nprobability', fontsize=14)

plt.savefig('./endog_stickplot.pdf', facecolor='white', bbox_inches='tight')

#%%
