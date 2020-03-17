#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.viz
import vdj.io
vdj.viz.plotting_style()

# %%
# Load the NHST data 
nhst = pd.read_csv('../../data/looping_frequency_p_values.csv')
dwell_time_nhst = pd.read_csv('../../data/dwell_time_p_values.csv')
nhst.rename(columns={'p_value':'floop_p'}, inplace=True)
nhst['dwell_p'] = dwell_time_nhst['p_value']

# Add positional labels to the mutants
nhst['position'] = np.arange(0, len(loop_freq_nhst), 1)
nhst.loc[nhst['floop_p'] <= 0.05, 'floop_color'] = 'dodgerblue'
nhst.loc[nhst['floop_p'] > 0.05, 'floop_color'] = '#3c3c3c'
nhst.loc[nhst['dwell_p'] <= 0.05, 'dwell_color'] = 'dodgerblue'
nhst.loc[nhst['dwell_p'] > 0.05, 'dwell_color'] = '#3c3c3c'

# %%
fig, ax = plt.subplots(1, 2, figsize=(5, 8.5), sharey=True)
for p in nhst['position'].values:
    if p%2 != 0:
        ax[0].hlines(p, -9E-6, 1.1, color='white', lw=10, alpha=0.5, zorder=1)
        ax[1].hlines(p, -9E-6, 1.1, color='white', lw=10, alpha=0.5, zorder=1)

for a in ax:
    a.set_xscale('symlog', linthreshx=1E-4)
    a.set_yticks(loop_freq_nhst['position'])
    a.set_yticklabels(loop_freq_nhst['mutant'])
    a.yaxis.set_tick_params(labelsize=6)
    a.set_ylim([-0.5, len(loop_freq_nhst) + 0.25])
    a.set_xlabel('$p$-value')
    a.set_xlim([-9E-6 ,1.2])

# Add titles
ax[0].set_title('Looping Frequency')
ax[1].set_title('Median Dwell Time')
for g, d in nhst.groupby(['mutant']):
    ax[0].plot(d['floop_p'], d['position'], 'o',
          markerfacecolor='white', color=d['floop_color'].values[0])
    if str(d['dwell_p'].values[0]) != 'nan':
        ax[1].plot(d['dwell_p'], d['position'], 'o',
          markerfacecolor='white', color=d['dwell_color'].values[0])
plt.savefig('../../figures/SiFigX_NHST.pdf', bbox_inches='tight')



# %%
