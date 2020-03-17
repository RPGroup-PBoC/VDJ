#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.viz
import vdj.io
vdj.viz.plotting_style()

# %%
# Looping Frequency NHST
floop_nhst = pd.read_csv('../../data/looping_frequency_p_values.csv')
floop_nhst.rename(columns={'p_value':'floop_p'}, inplace=True)

# Dwell Time NHST
dwell_time_nhst = pd.read_csv('../../data/dwell_time_p_values.csv')
dwell_time_nhst.rename(columns={'p_value':'dwell_p'}, inplace=True)

# Pcut NHST
pcut_nhst = pd.read_csv('../../data/cutting_probability_p_values.csv')
pcut_nhst.rename(columns={'p_value':'pcut_p'}, inplace=True)

nhst = floop_nhst.merge(dwell_time_nhst, on='mutant')
nhst = nhst.merge(pcut_nhst, on='mutant')

#%%
# Filter the mutants
muts = []
for g, d in nhst.groupby(['mutant']):
    parsed = vdj.io.mutation_parser(g)
    if ('Spac' in g) | ('Non' in g) | ('Hept' in g):
        if parsed['n_muts'] == 1:
            muts.append(g)
    elif 'Cod' not in g:
        muts.append(g)
nhst = nhst[nhst['mutant'].isin(muts)]
#%%

# Add positional labels to the mutants
nhst['position'] = np.arange(0, len(nhst), 1)
nhst.loc[nhst['floop_p'] <= 0.05, 'floop_color'] = 'dodgerblue'
nhst.loc[nhst['floop_p'] > 0.05, 'floop_color'] = '#3c3c3c'
nhst.loc[nhst['dwell_p'] <= 0.05, 'dwell_color'] = 'dodgerblue'
nhst.loc[nhst['dwell_p'] > 0.05, 'dwell_color'] = '#3c3c3c'
nhst.loc[nhst['pcut_p'] <= 0.05, 'pcut_color'] = 'dodgerblue'
nhst.loc[nhst['pcut_p'] > 0.05, 'pcut_color'] = '#3c3c3c'


# %%
fig, ax = plt.subplots(1, 3, figsize=(6, 8.5), sharey=True)
for p in nhst['position'].values:
    if p%2 != 0:
        ax[0].hlines(p, -9E-6, 1.1, color='white', lw=10, alpha=0.5, zorder=1)
        ax[1].hlines(p, -9E-6, 1.1, color='white', lw=10, alpha=0.5, zorder=1)
        ax[2].hlines(p, -9E-6, 1.1, color='white', lw=10, alpha=0.5, zorder=1)

for a in ax:
    a.set_xscale('symlog', linthreshx=1E-4)
    a.set_yticks(nhst['position'])
    a.set_yticklabels(nhst['mutant'])
    a.yaxis.set_tick_params(labelsize=6)
    a.xaxis.set_tick_params(labelsize=6)
    a.set_ylim([-0.5, len(nhst) + 0.25])
    a.set_xlabel('$p$-value')
    a.set_xlim([-9E-6 ,1.2])

# Add titles
ax[0].set_title('Looping Frequency')
ax[1].set_title('Median Dwell Time')
ax[2].set_title('Cutting Probability')
for g, d in nhst.groupby(['mutant']):
    ax[0].plot(d['floop_p'], d['position'], 'o',
          markerfacecolor='white', color=d['floop_color'].values[0])
    ax[1].plot(d['dwell_p'], d['position'], 'o',
          markerfacecolor='white', color=d['dwell_color'].values[0])
    ax[2].plot(d['pcut_p'], d['position'], 'o',
          markerfacecolor='white', color=d['pcut_color'].values[0])

plt.savefig('../../figures/SiFigX_NHST.pdf', bbox_inches='tight')



# %%
