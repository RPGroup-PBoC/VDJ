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
map = {"DFL16.1-3'":1, "DFL16.1-5'":2, 'V1-135':3, 'V9-120':4,
       'V10-96':5, 'V19-93':6, 'V4-57-1 (ref)':7, 'V4-55':8,
       'V5-43':9, 'V8-18':10, 'V6-17':11, 'V6-15':12}
muts = list(map.keys())
#%%
post_colors = {'12HeptC3G' : '#E10C00', 
                '12HeptC3T' : '#BF3030', 
                '12SpacC4G' : '#A24F59', 
                '12NonA4C' : '#7D778E', 
                'V4-57-1 (ref)' :  'slategrey', #599DC1', 
                '12SpacG10T' : '#679CE8'} #'#38C2F2'}

posterior_offset = np.linspace(0, 0.025 * (len(muts) - 1), len(muts))
offset_dict = {m:p for m, p in zip(muts, posterior_offset)}

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
for mut, d in endo_cuts.groupby(['mutant']):
    ax.fill_between(d['probability'], 0 + offset_dict[mut], 
                    d['posterior'] + offset_dict[mut],
                    zorder=len(muts) - map[mut], lw=0)
    ax.plot(d['probability'], d['posterior'] + offset_dict[mut],
            zorder=len(muts) - map[mut])
    ax.text(1.0, offset_dict[mut] - 0.005, mut, 
            zorder=20, ha='right')
ax.set_xlim([0.0, 1.0])
ax.set_xlabel('$p_{cut}$', fontsize=14)
#%%
