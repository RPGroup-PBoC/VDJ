#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

# Load the data with long-form looping confidence intervals 
# and restrict to relevant sets.
data = pd.read_csv('../../data/compiled_loop_freq_bs.csv')
data = data[(data['hmgb1']==80)]

camg_data = pd.DataFrame()

for mut in data[data['salt']=='Ca']['mutant'].unique():
    camg_data = camg_data.append(data[data['mutant']==mut], ignore_index=True)
camg_data = camg_data.replace(to_replace='WT12rss', value='reference')
#%%
# Provide three different plots for calcium-magnesium looping
# frequency plots
muts = {'reference':0, '12HeptA4T':1, '12SpacG11T':2}
df_loops = camg_data[camg_data['percentile']==95.0].copy()

fig, ax = plt.subplots(3, 1, figsize=(6,6))
for mut in muts:
    ca_loops = df_loops[(df_loops['mutant']==mut) & (df_loops['salt']=='Ca')]['loops_per_bead'].values[0]
    mg_loops = df_loops[(df_loops['mutant']==mut) & (df_loops['salt']=='Mg')]['loops_per_bead'].values[0]
    ax[muts[mut]].set_xlim([0, 1.0])
    ax[muts[mut]].set_ylim([0, 0.1])
    ax[muts[mut]].set_yticklabels([])

    if muts[mut]!=2:
        ax[muts[mut]].set_xticklabels([])

    ax[muts[mut]].scatter(ca_loops, 0.06, color='green', label=r'Ca$^{2+}$')
    ax[muts[mut]].scatter(mg_loops, 0.03, color='rebeccapurple', label=r'Mg$^{2+}$')
    
    for g,d in camg_data[camg_data['mutant']==mut].groupby(['percentile','salt']):
        if g[1]=='Ca':
            rect = patches.Rectangle([d['low'].values[0], 0.048], 
                                    (d['high'].values[0] - d['low'].values[0]),
                                    0.024, color='green', alpha=0.15)
            ax[muts[mut]].add_patch(rect)
        elif g[1]=='Mg':
            rect = patches.Rectangle([d['low'].values[0], 0.018],
                                    (d['high'].values[0] - d['low'].values[0]),
                                    0.024, color='rebeccapurple', alpha=0.15)
            ax[muts[mut]].add_patch(rect)
    ax[muts[mut]].text(0.98, 0.085, mut, fontsize=12, ha='right')
ax[2].set_xlabel('loop frequency', fontsize=14)

# Draw rectangles as legend
confidence_intervals = np.sort(camg_data['percentile'].unique())
widths = np.linspace(0, 0.02 * (len(confidence_intervals) - 1), len(confidence_intervals))
ci_dict = {str(ci):w for ci,w in zip(confidence_intervals, widths)}
for ci in ci_dict:
    ci_rect_ca = patches.Rectangle((0.2, 0.105), ci_dict[ci], 0.024, 
                                    color='green', alpha=0.15, clip_on=False)
    ci_rect_mg = patches.Rectangle((0.65, 0.105), ci_dict[ci], 0.024, 
                                    color='rebeccapurple', alpha=0.15, clip_on=False)
    ax[0].add_patch(ci_rect_ca)
    ax[0].add_patch(ci_rect_mg)

#%%


#%%
