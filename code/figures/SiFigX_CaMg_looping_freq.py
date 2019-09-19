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
    ax[muts[mut]].set_ylim([0.01, 0.08])
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
    ax[muts[mut]].text(0.98, 0.015, mut, fontsize=14, ha='right')
ax[2].set_xlabel('loop frequency', fontsize=14)

# Draw rectangles as legend
confidence_intervals = np.sort(camg_data['percentile'].unique())
widths = np.linspace(0.05, 0.06 * (len(confidence_intervals)), len(confidence_intervals))
ci_dict = {str(int(ci)):w for ci,w in zip(confidence_intervals, widths)}
for ci in ci_dict:
    ci_rect_ca = patches.Rectangle((0.10, 0.095), ci_dict[ci], 0.024, 
                                    color='green', alpha=0.15, clip_on=False)
    ci_rect_mg = patches.Rectangle((0.60, 0.095), ci_dict[ci], 0.024, 
                                    color='rebeccapurple', alpha=0.15, clip_on=False)
    ax[0].add_patch(ci_rect_ca)
    ax[0].add_patch(ci_rect_mg)
    ax[0].text(0.08 + ci_dict[ci], 0.085, ci + '%', ha='center')
    ax[0].text(0.58 + ci_dict[ci], 0.085, ci + '%', ha='center')
ax[0].text(0.25, 0.125, r'Ca$^{2+}$', fontsize=14)
ax[0].text(0.75, 0.125, r'Mg$^{2+}$', fontsize=14)
plt.tight_layout()
plt.savefig('./SiFigX_CaMg_looping_freq.pdf', facecolor='white', bbox_inches='tight')
#%%


#%%
