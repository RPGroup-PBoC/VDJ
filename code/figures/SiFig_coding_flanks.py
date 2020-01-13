"""
Comparison of Coding Flank Effects and Critical SpacerC1A Mutant
--------------------------------------------------------------------------------
Author: Soichi Hirokawa
Last Modified: January 7, 2020
License: MIT

Description
--------------------------------------------------------------------------------
This script generates SI Figures 4 and 5 which shows the effect of variable
coding flanks on the looping frequency, dwell time, and cutting probability
on the reference V4-57-1 RSS as well as effects from a C to A mutation at the
first positon of the spacer.

Notes
--------------------------------------------------------------------------------
This script is designed to be executed from the `code/figures` directory and
loads the proper CSV files from a relative path.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

# Upload data on coding flank-relevant RSSs
cf_muts = ['WT12rss', '12CodC6A', '12SpacC1A', 'V4-55']

loop = pd.read_csv('../../data/compiled_looping_frequency_bootstrap.csv',
                   comment='#')
loop = loop[(loop['mutant'].isin(cf_muts)) &
         (loop['hmgb1']==80) & (loop['salt']=='Mg')]
loop = loop.replace(to_replace='WT12rss', value='V4-57-1 (ref)')

dwell = pd.read_csv('../../data/compiled_dwell_times.csv', comment='#')
dwell = dwell[(dwell['mutant'].isin(cf_muts)) & (dwell['hmgb1']==80) & (dwell['salt']=='Mg')]
dwell = dwell.replace(to_replace='WT12rss', value='V4-57-1 (ref)')
dwell_cut = dwell[dwell['cut']==1].copy()
dwell_unloop = dwell[dwell['cut']==0].copy()

cut_posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv',
                        comment='#')
cut_posts = cut_posts[(cut_posts['mutant'].isin(cf_muts)) & (cut_posts['hmgb1']==80) & (cut_posts['salt']=='Mg')]
cut_posts = cut_posts.replace(to_replace='WT12rss', value='V4-57-1 (ref)')

# Create dwell ECDFs
dfs = []
for source in [dwell, dwell_cut, dwell_unloop]:
    dist_df = []
    for g,d in source.groupby('mutant'):
        x,y = np.sort(d['dwell_time_min'].values), np.arange(0, len(d), 1) / len(d)
        y[-1] = 1
        _df = pd.DataFrame()
        _df['x'] = x
        _df['y'] = y
        _df['mutant'] = g
        dist_df.append(_df)
    dwell_dist = pd.concat(dist_df)
    dfs.append(dwell_dist)
dwell_dist, cut_dist, unloop_dist = dfs

#%%
fig = plt.figure()
ax = plt.gca()
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.set_facecolor('white')

gs_ci = gridspec.GridSpec(1, 3)
ax_conf_int = fig.add_subplot(gs_ci[:, :])
gs_ci.tight_layout(fig, rect=[0, 1.5, 1.5, 1.75])
ax_conf_int.set_xlim([0.15, 0.78])
ax_conf_int.set_xticklabels([])
ax_conf_int.set_yticklabels([])
ax_conf_int.set_facecolor('white')

gs = gridspec.GridSpec(1, 3, hspace=0.2, wspace=0.6)
ax_loop = fig.add_subplot(gs[0])
ax_loop.set_xticklabels([])
ax_loop.set_xlim([0.125, 0.625])
ax_loop.set_ylim([0, 0.3])
ax_loop.set_yticklabels([0, '', 0.1, '', 0.2, '', 0.3])
ax_loop.set_ylabel('looping frequency', fontsize=16)

ax_posts = fig.add_subplot(gs[1:])
ax_posts.set_xlabel('cutting probability', fontsize=16)
ax_posts.set_ylabel('posterior distribution', fontsize=16)

gs.tight_layout(fig, rect=[0, 0.75, 1.5, 1.55])

gs2 = gridspec.GridSpec(1, 3)
ax_dwell_cut = fig.add_subplot(gs2[0])
ax_dwell_cut.set_xscale('log')
ax_dwell_cut.text(14, 0.1, 'cutting\nonly',
                horizontalalignment='center', fontsize=14)
ax_dwell_cut.set_xlabel('dwell time [min]', fontsize=16)
ax_dwell_cut.set_ylabel('ECDF', fontsize=16)

ax_dwell_unloop = fig.add_subplot(gs2[1])
ax_dwell_unloop.set_xscale('log')
ax_dwell_unloop.text(14, 0.1, 'unlooping\nonly',
                horizontalalignment='center', fontsize=14)
ax_dwell_unloop.set_xlabel('dwell time [min]', fontsize=16)
ax_dwell_unloop.set_yticklabels([])

ax_dwell_all = fig.add_subplot(gs2[2])
ax_dwell_all.set_xscale('log')
ax_dwell_all.text(14, 0.1, 'cutting and\nunlooping',
                horizontalalignment='center', fontsize=14)
ax_dwell_all.set_xlabel('dwell time [min]', fontsize=16)
ax_dwell_all.set_yticklabels([])

gs2.tight_layout(fig, rect=[0, 0, 1.5, 0.8])

seqs = {'V4-57-1 (ref)':0.25, '12CodC6A':0.5}
colors = {'V4-57-1 (ref)':'slategrey', '12CodC6A':'dodgerblue'}
perc_widths = {p:w for p, w in zip(np.sort(loop['percentile'].unique()), np.arange(0.03, 0.21, 0.03))}

for seq in seqs:
    ax_loop.scatter(seqs[seq], loop[loop['mutant']==seq]['loops_per_bead'].values[0], 
                    s=150, marker='^', color='white', edgecolor=colors[seq], zorder=10)
    for g,d in loop[loop['mutant']==seq].groupby('percentile'):
        rect = patches.Rectangle([seqs[seq] - 0.075, d['low']], 0.15,
                                d['high'].values[0] - d['low'].values[0],
                                color=colors[seq], alpha=0.2)
        ax_loop.add_patch(rect)
        rect2 = patches.Rectangle([seqs[seq], 0.1], perc_widths[g], 0.5, alpha=0.2, color=colors[seq])
        ax_conf_int.add_patch(rect2)
        ax_conf_int.text(seqs[seq] + perc_widths[g] - 0.015, 0.05, str(int(g)) + '%',
                        ha='center', va='top')

    ax_dwell_cut.step(cut_dist[cut_dist['mutant']==seq]['x'], 
                        cut_dist[cut_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_dwell_unloop.step(unloop_dist[unloop_dist['mutant']==seq]['x'], 
                        unloop_dist[unloop_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_dwell_all.step(dwell_dist[dwell_dist['mutant']==seq]['x'], 
                        dwell_dist[dwell_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_posts.plot(cut_posts[cut_posts['mutant']==seq]['probability'],
                    cut_posts[cut_posts['mutant']==seq]['posterior'],
                    lw=2, color=colors[seq])
    ax_posts.fill_between(cut_posts[cut_posts['mutant']==seq]['probability'],
                            0, cut_posts[cut_posts['mutant']==seq]['posterior'],
                            color=colors[seq], alpha=0.4)
    
    ax_conf_int.text(seqs[seq] + 0.09, 0.65, seq, ha='center', va='bottom', fontsize=12)
plt.savefig('../../figures/SiFig_coding_flank_reference.pdf', bbox_inches='tight', facecolor='white')
#%%
seqs = {'V4-55':0.25, '12SpacC1A':0.5}
colors = {'V4-55':'slategrey', '12SpacC1A':'dodgerblue'}

fig = plt.figure()
ax = plt.gca()
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ax.set_facecolor('white')

gs_ci = gridspec.GridSpec(1, 3)
ax_conf_int = fig.add_subplot(gs_ci[:, :])
gs_ci.tight_layout(fig, rect=[0, 1.5, 1.5, 1.75])
ax_conf_int.set_xlim([0.15, 0.78])
ax_conf_int.set_xticklabels([])
ax_conf_int.set_yticklabels([])
ax_conf_int.set_facecolor('white')

gs = gridspec.GridSpec(1, 3, hspace=0.2, wspace=0.6)
ax_loop = fig.add_subplot(gs[0])
ax_loop.set_xticklabels([])
ax_loop.set_xlim([0.125, 0.625])
ax_loop.set_ylim([0, 0.3])
ax_loop.set_yticklabels([0, '', 0.1, '', 0.2, '', 0.3])
ax_loop.set_ylabel('looping frequency', fontsize=16)

ax_posts = fig.add_subplot(gs[1:])
ax_posts.set_xlabel('cutting probability', fontsize=16)
ax_posts.set_ylabel('posterior distribution', fontsize=16)

gs.tight_layout(fig, rect=[0, 0.75, 1.5, 1.55])

gs2 = gridspec.GridSpec(1, 3)
ax_dwell_cut = fig.add_subplot(gs2[0])
ax_dwell_cut.text(12, 0.1, 'cutting\nonly', 
                horizontalalignment='center', fontsize=14)
ax_dwell_cut.set_xscale('log')
ax_dwell_cut.set_xlabel('dwell time [min]', fontsize=16)
ax_dwell_cut.set_ylabel('ECDF', fontsize=16)

ax_dwell_unloop = fig.add_subplot(gs2[1])
ax_dwell_unloop.text(12, 0.1, 'unlooping\nonly', 
                horizontalalignment='center', fontsize=14)
ax_dwell_unloop.set_xscale('log')
ax_dwell_unloop.set_xlabel('dwell time [min]', fontsize=16)
ax_dwell_unloop.set_yticklabels([])

ax_dwell_all = fig.add_subplot(gs2[2])
ax_dwell_all.text(10, 0.1, 'cutting and\nunlooping',
                horizontalalignment='center', fontsize=14)
ax_dwell_all.set_xscale('log')
ax_dwell_all.set_xlabel('dwell time [min]', fontsize=16)
ax_dwell_all.set_yticklabels([])

gs2.tight_layout(fig, rect=[0, 0, 1.5, 0.8])

perc_widths = {p:w for p, w in zip(np.sort(loop['percentile'].unique()), np.arange(0.03, 0.21, 0.03))}

for seq in seqs:
    ax_loop.scatter(seqs[seq], loop[loop['mutant']==seq]['loops_per_bead'].values[0], 
                    s=150, marker='^', color='white', edgecolor=colors[seq], zorder=10)
    for g,d in loop[loop['mutant']==seq].groupby('percentile'):
        rect = patches.Rectangle([seqs[seq] - 0.075, d['low']], 0.15,
                                d['high'].values[0] - d['low'].values[0],
                                color=colors[seq], alpha=0.2)
        ax_loop.add_patch(rect)
        rect2 = patches.Rectangle([seqs[seq], 0.1], perc_widths[g], 0.5, alpha=0.2, color=colors[seq])
        ax_conf_int.add_patch(rect2)
        ax_conf_int.text(seqs[seq] + perc_widths[g] - 0.015, 0.05, str(int(g)) + '%',
                        ha='center', va='top')

    ax_dwell_cut.step(cut_dist[cut_dist['mutant']==seq]['x'], 
                        cut_dist[cut_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_dwell_unloop.step(unloop_dist[unloop_dist['mutant']==seq]['x'], 
                        unloop_dist[unloop_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_dwell_all.step(dwell_dist[dwell_dist['mutant']==seq]['x'], 
                        dwell_dist[dwell_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_posts.plot(cut_posts[cut_posts['mutant']==seq]['probability'],
                    cut_posts[cut_posts['mutant']==seq]['posterior'],
                    lw=2, color=colors[seq])
    ax_posts.fill_between(cut_posts[cut_posts['mutant']==seq]['probability'],
                            0, cut_posts[cut_posts['mutant']==seq]['posterior'],
                            color=colors[seq], alpha=0.4)
    
    ax_conf_int.text(seqs[seq] + 0.09, 0.65, seq, ha='center', va='bottom', fontsize=12)
plt.savefig('../../figures/SiFig_spacerC1A_endogenous_cmparison.pdf', 
            bbox_inches='tight', facecolor='white')

# %%
