# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import imp
imp.reload(vdj.io)
import vdj.io
import vdj.stats
bokeh.io.output_notebook()

#%% 
# Load the sampling information
cut_stats = pd.read_csv('../../data/pooled_cutting_probability_summary.csv')

# Get the reference sequence
seqs = vdj.io.endogenous_seqs()
nt_idx = vdj.io.nucleotide_idx()
ref_seq  = seqs['reference'][0]
ref_int = seqs['reference'][1]

# Define the colors
colors = {'A': '#ee6352', 'T':'#08b2e3', 'C':'#57a773', 'G':'#484d6d'}

# Define the jitter
jitter = {'A': -0.1, 'T':-0.05, 'C':0.05, 'G':0.1}
# %%
# isolate the point mutants

# Define the nucleotide position
for m in cut_stats['mutant'].unique():
    # Parse the mutation 
    mut = vdj.io.mutation_parser(m)
    ind = np.argmax(mut['seq_idx'] != ref_int)
    idx = nt_idx[mut['seq'][ind]]
    base  = nt_idx[mut['seq_idx'][ind]]

    # Add the identity to the dataframe 
    cut_stats.loc[cut_stats['mutant']==m, 'n_muts'] = mut['n_muts']
    cut_stats.loc[cut_stats['mutant']==m, 'mut_idx'] = idx
    cut_stats.loc[cut_stats['mutant']==m, 'mut_base'] = base
    cut_stats.loc[cut_stats['mutant']==m, 'position'] = ind + jitter[base]

    # Add the plotting information
    cut_stats.loc[cut_stats['mutant']==m, 'color'] = colors[base]

point_muts = cut_stats[cut_stats['n_muts']==1]

wt = cut_stats[cut_stats['mutant']=='WT12rss']
# %%
ax = bokeh.plotting.figure(height=300, width=900, y_axis_label='cleavage probability',
                          x_axis_label='reference sequence', 
                          y_range=[-0.05, 1.05], x_range=[0, 29])
ax.quad(left=0, right=29, bottom=wt['hpd_min'], top=wt['hpd_max'], color='navy', alpha=0.2, legend='WT')
ax.segment(x0='position', x1='position', y0='hpd_min', y1='hpd_max', line_width=1, 
           color='color', source=point_muts)
ax.circle(x='position', y='median', line_color='color', fill_color='white', source=point_muts, legend='mut_base') 
ax.xaxis.ticker = np.arange(1, 29)
ax.xaxis.major_label_overrides = {str(i):a for i, a in zip(np.arange(1, 29), list(ref_seq))} 
ax.xgrid.grid_line_width = 0

x_bands = np.arange(1, 29, 2)
ax.quad(left=x_bands-0.45, right=x_bands+0.45, bottom=-0.05, top=1.05, color='gray', alpha=0.1)

ax.legend.location = 'top_right'
ax.legend.label_text_font_size = '8pt'
ax.legend.orientation = 'horizontal'
bokeh.io.show(ax)
bokeh.io.save(ax, filename='cutting_prob_pooled.html')
#%%

np.arange(1, 29, 2)




#%%
