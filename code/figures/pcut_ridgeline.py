# -*-coding: utf-8 -*-
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import vdj.io
import vdj.viz


# Load the posteriors
posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv')
points = posts[(posts['n_muts']==1) & (posts['salt']=='Mg') & (posts['hmgb1']==80)]

# Lod the reference identities
ref = vdj.io.endogenous_seqs()['WT12rss']
ref_seq = ref[0]
ref_idx = ref[1]

for m in points['mutant'].unique():
    seq = vdj.io.mutation_parser(m)
    loc = np.argmax(ref_idx != seq['seq_idx'])
    base = seq['seq'][loc]
    points.loc[points['mutant']==m, 'pos'] = loc
    points.loc[points['mutant']==m, 'base'] = base


# %%
# Set up the figure canvas
fig, ax = plt.subplots(28, 1, sharex=True, sharey=True, figsize=(3.42, 9))
for a in ax:
    a.set_yticks([])
    a.set_facecolor('none')
    a.hlines(0, 0, 1, 'k', lw=1)

# # Plot the wild-type posterior
# wt_post = posts[posts['mutant']=='WT12rss']
# for i in range(27):
#     ax[i].plot(wt_post['probability'], wt_post['posterior']/wt_post['posterior'].max(), lw=0.75, color='grey', alpha=0.75)
#     ax[i].fill_between(wt_post['probability'], wt_post['posterior']/wt_post['posterior'].max(), color='grey', alpha=0.1)
    
colors = {'A':'tomato', 'T':'dodgerblue', 'C':'grey', 'G':'purple'}
for g, d in points.groupby(['pos', 'base']):
    ax[int(g[0])].plot(d['probability'], d['posterior']/d['posterior'].max(), lw=0.75, color=colors[g[1]])
    ax[int(g[0])].fill_between(d['probability'], d['posterior']/d['posterior'].max(), lw=0.75, color=colors[g[1]],
                        alpha=0.1)

#%%
