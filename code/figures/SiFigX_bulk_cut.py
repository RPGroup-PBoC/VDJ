#%%
"""
Generates SI Figure for bulk DNA tether loss with data on looping fraction
and cutting probability.
"""
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import vdj.viz
import vdj.io
vdj.viz.plotting_style()

def mutation_value(mutant):
        base_val = {'Hept':0, 'Spac':30, 'Non':80}
        endog_val = {'DFL1613':1, 'DFL161':2, 'V1-135':3, 'V9-120':4,
                        'V10-96':5, 'V19-93':6, 'WT12rss':7, 'V4-55':8,
                        'V5-43':9, 'V8-18':10, 'V6-17':11, 'V6-15':12}
        nuc_val = {'A':0, 'C':1, 'G':2, 'T':3}
        if mutant[:2]!='12':
                # random large number for endogenous RSSs
                num = 200 + endog_val[mutant]
        elif 'Non' in mutant:
                num = base_val[mutant[2:5]] + 4*int(mutant[6:-1]) + nuc_val[mutant[-1]]
        else:
                num = base_val[mutant[2:6]] + 4*int(mutant[7:-1]) + nuc_val[mutant[-1]]
        return num
#%%
# Load the data with long-form looping events and restrict to relevant sets.
data = pd.read_csv('../../data/compiled_loop_freq_bs.csv',       
                   comment='#')
counts = data[(data['salt']=='Mg') & (data['hmgb1']==80) & 
                (data['percentile']==95.0) & (data['mutant']!='12CodC6A')]
counts = counts.drop(columns=['n_beads', 'n_loops', 'percentile'])

# Load all cutting probability estimates taking Gaussian approximation.
cut_data = pd.read_csv('../../data/pooled_cutting_probability.csv', comment='#')
cut_data = cut_data[(cut_data['hmgb1'] == 80) & (cut_data['salt']=='Mg') & 
                (cut_data['mutant']!='12CodC6A')]
cut_data.rename(columns={'mode':'cut_mode', 'std':'cut_std'},
                inplace=True)
cut_data = cut_data.drop(columns=['n_beads', 'n_cuts', 'n_muts'])

# Load all bead cleavage probability estimates taking Gaussian approximation.
bead_data = pd.read_csv('../../data/pooled_bead_cutting_probability.csv', comment='#')
bead_data = bead_data[(bead_data['hmgb1'] == 80) & (bead_data['salt']=='Mg') & 
                        (bead_data['mutant']!='12CodC6A')]
bead_data.rename(columns={'mode':'bead_mode', 'std':'bead_std'},
                inplace=True)

df = pd.merge(bead_data, counts, on=['mutant', 'hmgb1', 'salt'])
df = pd.merge(df, cut_data, on=['mutant', 'hmgb1', 'salt'])
df = df.drop(columns=['Unnamed: 0'])

# Filter the mutants
muts = []
for g, d in df.groupby(['mutant']):
    parsed = vdj.io.mutation_parser(g)
    if ('Spac' in g) | ('Non' in g) | ('Hept' in g):
        if parsed['n_muts'] == 1:
            muts.append(g)
    elif 'Cod' not in g:
        muts.append(g)
df = df[df['mutant'].isin(muts)]
df = df[df['n_muts']!=2.0]

df['ordering'] = 0
for mut in df['mutant']:
        df.loc[(df['mutant']==mut), 'ordering'] = mutation_value(mut)

df = df.sort_values('ordering', ignore_index=True)
df = df.replace({'DFL161':"DFL16.1-5'", 'DFL1613':"DFL16.1-3'",
                'WT12rss':'V4-57-1 (ref)'})
wt = df[df['mutant']=='V4-57-1 (ref)']

mutants = df[df['mutant']!='V4-57-1 (ref)']
mutants['position'] = np.arange(0, len(mutants), 1)
#%%
fig, ax = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
for p in mutants['position'].values:
    if p%2 != 0:
        ax[0].vlines(p, 0, 1.0, color='white', lw=12, alpha=0.5, zorder=1)
        ax[1].vlines(p, 0, 1.0, color='white', lw=12, alpha=0.5, zorder=1)
        ax[2].vlines(p, 0, 1.0, color='white', lw=12, alpha=0.5, zorder=1)

xlabel = ['loops per bead', 'cutting probability', 'bead cut fraction']

for n in range(len(ax)):
        ax[n].set_xticks(mutants['position'])
        ax[n].set_xticklabels(mutants['mutant'])
        ax[n].yaxis.set_tick_params(labelsize=10)
        ax[n].xaxis.set_tick_params(labelsize=10, rotation=90)
        ax[n].set_xlim([-0.5, len(mutants)-0.5])
        ax[n].set_ylabel('%s' %xlabel[n], fontsize=16)
ax[0].set_ylim([0,0.65])
ax[1].set_ylim([0,1.0])
ax[2].set_ylim([0,0.3])

# Add plots with staggering
metric = ['loops_per_bead', 'cut_mode', 'bead_mode']
lower = ['low', 'cut_std', 'bead_std']
upper = ['high', 'cut_std', 'bead_std']
label = ['loops per bead', 'cutting probability', 'bead cleavage fraction']
colors = ['tomato', 'dodgerblue', 'rebeccapurple']

for g, d in mutants.groupby(['mutant']):
        for n in range(len(ax)):
                ax[n].plot(d['position'], d[metric[n]], 'o',
                        markerfacecolor='white', color=colors[n])
                if n==0:
                        ax[n].errorbar(d['position'], d[metric[n]], 
                                yerr=(d[metric[n]]-d[lower[n]],d[upper[n]]-d[metric[n]]),
                                lw=2, color=colors[n])
                else:
                        ax[n].errorbar(d['position'], d[metric[n]], 
                                yerr=(d[lower[n]],d[upper[n]]),
                                lw=2, color=colors[n])
for n in range(len(ax)):
        ax[n].axhline(wt[metric[n]].values[0], -0.5, len(mutants)-0.5,
                        color='k', ls='--', zorder=1)
        if n==0:
                ax[n].axhspan(wt[lower[n]].values[0], wt[upper[n]].values[0],
                        color='k', alpha=0.5, zorder=1)
        else:
                ax[n].axhspan(wt[metric[n]].values[0]-wt[lower[n]].values[0],
                                wt[upper[n]].values[0]+wt[metric[n]].values[0],
                                color='k', alpha=0.5, zorder=1)

plt.savefig('./SiFigX_bead_cut.pdf', bbox_inches='tight', facecolor='white')

# %%
