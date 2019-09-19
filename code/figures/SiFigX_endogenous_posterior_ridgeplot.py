#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

# Load the precomputed posterior distributioons
cut_posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv')
cut_posts = cut_posts[(cut_posts['hmgb1']==80) & (cut_posts['salt']=='Mg')]

#%%
df = []
for g, d in cut_posts.groupby(['mutant']):
    if ('Hept' not in g) & ('Spac' not in g) & ('Non' not in g):
        df.append(d) 

endo_cuts = pd.concat(df)
mut_names = ['DFL161', 'DFL1613', 'WT12rss']
new_names = ["DFL16.1-5'", "DFL16.1-3'", 'V4-57-1 (ref)']
endo_cuts = endo_cuts.replace(to_replace=mut_names, value=new_names)
muts = endo_cuts['mutant'].unique()

map = {m:i for m, i in zip(muts, np.arange(1, 13))}
map = {'V6-15':12, 'V6-17':11, 'V8-18':10, 'V5-43':9,
       'V4-55':8, 'V4-57-1 (ref)':7, 'V19-93':6, 'V10-96':5,
       'V9-120':4, 'V1-135':3, "DFL16.1-5'":2, "DFL16.1-3'":1}
muts = list(map.keys())
#%%
post_colors = {"DFL16.1-3'" : '#2ecc71', 
                "DFL16.1-5'" : '#52be80', 
                'V1-135' : '#21618c', 
                'V9-120' : '#2874a6', 
                'V10-96' : '#3498db',
                'V19-93' : '#5dade2',
                'V4-57-1 (ref)' :  'slategrey', 
                'V4-55' : '#ec7063',
                'V5-43' : '#e74c3c',
                'V8-18' : '#cb4335',
                'V6-17' : '#b03a2e',
                'V6-15' : '#943126'}

posterior_offset = np.linspace(0, 0.025 * (len(muts) - 1), len(muts))
offset_dict = {m:p for m, p in zip(muts, posterior_offset)}

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
for mut, d in endo_cuts.groupby(['mutant']):
    ax.fill_between(d['probability'], 0 + offset_dict[mut], 
                    d['posterior'] + offset_dict[mut],
                    zorder=len(muts) - map[mut], lw=0, alpha=0.4,
                    facecolor=post_colors[mut])
    ax.plot(d['probability'], d['posterior'] + offset_dict[mut],
            zorder=len(muts) - map[mut], color=post_colors[mut], alpha=0.6)
    ax.text(1.0, offset_dict[mut] - 0.011, mut, 
            zorder=20, ha='right', color=post_colors[mut])
ax.vlines(0.73, offset_dict["DFL16.1-3'"], offset_dict["DFL16.1-3'"] + 0.025,
        color='k', lw=2)
ax.hlines(offset_dict["DFL16.1-3'"], 0.7, 0.73, lw=2)
ax.hlines(offset_dict["DFL16.1-3'"] + 0.025, 0.7, 0.73, lw=2)
ax.text(0.74, offset_dict["DFL16.1-3'"] + 0.012, '$\propto$probability', fontsize=12)
ax.set_xlim([0.0, 1.0])
ax.set_yticklabels([])
ax.set_xlabel('$p_{cut}$', fontsize=14)
plt.savefig('./SiFigX_endogenous_posterior_ridgeplot.pdf',
            facecolor='white', bbox_inches='tight')
#%%
