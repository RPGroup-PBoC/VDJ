# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cor
import matplotlib.cm as cm
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

# Load the statistics and posterior data sets
stats = pd.read_csv('../../data/cutting_rate_analytic_summary.csv')
posteriors = pd.read_csv('../../data/cutting_rate_posteriors.csv')

# Isolate the wild-type stats
wt_stats = stats[stats['mutant']=='WT12rss']

# Keep only the point mutants
points = stats[(stats['n_muts'] == 1) & (stats['mutant'] != 'V4-55')].copy()
posts = posteriors[(posteriors['n_muts'] == 1) & (posteriors['mutant'] != 'V4-55')].copy()
# Load the reference sequence
seqs = vdj.io.endogenous_seqs()
ref_seq = seqs['WT12rss'][0]
ref_idx = seqs['WT12rss'][1]

# Load the nucleotide identities
nt_idx = vdj.io.nucleotide_idx()

wt_err_high = wt_stats['95_high'] / wt_stats['median']
wt_err_low = wt_stats['95_low'] / wt_stats['median']
wt_width = wt_err_high - wt_err_low

# Iterate through each point mutant and determine the location. 
for g, d in points.groupby(['mutant']):
    seq = vdj.io.mutation_parser(g)
    idx = seq['seq_idx']
    loc = np.argmax(ref_idx !=  idx)
    base = seq['seq'][loc]

    # Add the positioning back in
    points.loc[points['mutant']==g, 'position'] = loc
    points.loc[points['mutant']==g, 'mutation'] = nt_idx[base] 

    # Compute the size of the bubble based on 1/log10of the width
    size = 7 * (d['median'] / (d['width'] / 2)) + 11 
    points.loc[points['mutant']==g, 'rate_size'] = size

    # Compute the difference to the wild-type value
    diff = d['median'].values[0] / wt_stats['median'].values[0]
    inf_err_high = d['95_high'] / d['median'].values[0]
    inf_err_low = d['95_low'] / d['median'].values[0]
    inf_err = inf_err_high - inf_err_low
               
    # points.loc[points['mutant']==g, 'relative_diff'] = diff
    # points.loc[points['mutant']==g, 'diff_width'] = np.abs(width)
    # points.loc[points['mutant']==g, 'diff_size'] = 15 * (wt_width.values[0] / inf_err.values[0])

# %%
# Define the colors
rate_norm = cor.Normalize(vmin=0, vmax=90)
diff_norm = cor.Normalize(vmin=-2, vmax=2)
blues = cm.plasma
puor = cm.PuOr

# %%
# ##############################################################################
#  FIGURE INSTANTIATION
# ##############################################################################
# Instantiate the figure canvas
fig, ax = plt.subplots(2, 1, figsize=(7.5, 5))

# Format the axes
for i in range(len(ax)):
    if i!=1:
        ax[i].set_xticks(np.arange(1, 30))
        ax[i].set_xticklabels(list(ref_seq))
        ax[i].set_yticks(np.arange(1, 5))
        ax[i].set_yticklabels([nt_idx[k] for k in range(4)])
        ax[i].yaxis.grid(True)
        ax[i].set_ylim([0.5, 5])
        ax[i].set_xlim([0.9, 28.4])
        ax[i].set_facecolor('white')
        ax[i].yaxis.grid(color='grey', linewidth=0.25)

ax[0].set_xlabel('wild-type sequence', fontsize=12)
ax[0].set_ylabel('mutation', fontsize=12)
ax[1].set_xlim([0, 60])
ax[1].set_yticks([])
ax[1].set_ylabel('$\propto$ probability', fontsize=12)
ax[1].set_xlabel('cleavage waiting time [min]', fontsize=12)

# Fill in the blank positions with white x's
for a in [ax[0]]:
    for i in range(0, 28):
        for j in range(0, 4):
            if len(points[(points['position'] == i) & (points['mutation']==j)])==0:
                a.plot(i+1, j + 1, 'x', color='grey', ms=5, markeredgewidth=0.25,
                label='__nolegend__')

ax[0].plot([], [], marker='o', ms=10 * (100/50), label=' ', linewidth=0, color='k')
ax[0].plot([], [], marker='o', ms=10 * (100/90), label=' ', linewidth=0, color='k')
ax[0].plot([], [], marker='o', ms=10 * (0.1), label=' ', linewidth=0, color='k')
ax[0].legend(title=r'$-$ decreasing confidence $\rightarrow$', ncol=3, bbox_to_anchor=(0.95, 1.3)) 
    
# ##############################################################################
#  RATE VALUE BUBBLE PLOT
# ##############################################################################
wt_post = posteriors[posteriors['mutant']=='WT12rss']
ax[1].plot(wt_post['tau'], wt_post['posterior_pdf'], '--', lw=2, color='k', label='wild-type')



for i in range(len(points)):
    p = points.iloc[i]
    if p['mutant'] in ['12HeptG5A', '12SpacT9G', '12SpacG11T']:
        ax[0].plot(p['position'] + 1, p['mutation'] + 1, marker='s', markeredgecolor='tomato', color='None', ms=25)
        dist = posts[posts['mutant']==p['mutant']]
        ax[1].plot(dist['tau'], dist['posterior_pdf'], lw=2, label=p['mutant'], 
                   color=blues(rate_norm(p['median'])))
        ax[1].fill_between(dist['tau'], dist['posterior_pdf'], 
                   color=blues(rate_norm(p['median'])), alpha=0.35)

    ax[0].plot(p['position']+1, p['mutation'] + 1, '.', ms=p['rate_size'], color=blues(rate_norm(p['median'])),
               markeredgecolor='black', markeredgewidth=1)
    # ax[-1].plot(p['position']+1, p['mutation'] + 1, '.', ms=p['rate_size'], color=puor(diff_norm(p['relative_diff'] - 1)),
            #    markeredgecolor='black', markeredgewidth=1)

    ax[1].legend(loc='upper right', fontsize=12)

plt.tight_layout()
plt.savefig('./cutting_rates.pdf', bbox_inches='tight')
#%%


#%%
