#%%
import numpy as np
import pandas as pd
import vdj.io
import vdj.bayes
import vdj.viz
import bokeh.io
import bokeh.plotting
import bokeh.palettes
import bokeh.transform
import bokeh.models
import seaborn as sns
import tqdm
import glob
vdj.viz.plotting_style()
import imp
bokeh.io.output_notebook()
# %%
# Load the data sets
pcut_stats = pd.read_csv('../../data/cutting_probability_summary.csv')
ploop_stats = pd.read_csv('../../data/looping_probability_summary.csv')

# %%
ref_seq = vdj.io.endogenous_seqs()['reference'][1]

# Get the cutting probability of the reference sequence
for s in [pcut_stats, ploop_stats]:
    ref_prob = s[s['mutant']=='WT12rss']['median'].values[0]
    s['relative_prob'] = s['median'].values[0] - ref_prob

# %%
# Trim the data set to only the point muations
cut_muts = pcut_stats[pcut_stats['n_muts']<=1]
loop_muts = ploop_stats[ploop_stats['n_muts']<=1]

# Get the converter
nt_idx = vdj.io.nucleotide_idx()
ref_seq = vdj.io.endogenous_seqs()['reference'][1]
# Iterate through each mutant and identify the location of the changed value
for g, d in cut_muts.groupby('mutant'):
    # Get the sequence
    seq = vdj.io.mutation_parser(g)
    ints = seq['seq_idx']
    diff = (np.abs(ref_seq - ints) > 0)
    # Get the position and the change. 
    pos = np.argmax(diff)
    change = list(seq['seq'])[pos]

    # Add the position and change as x and y
    cut_muts.loc[cut_muts['mutant']==g, 'x'] = pos + 1
    cut_muts.loc[cut_muts['mutant']==g, 'y'] = nt_idx[change]

for g, d in loop_muts.groupby('mutant'):
    # Get the sequence
    seq = vdj.io.mutation_parser(g)
    ints = seq['seq_idx']
    diff = (np.abs(ref_seq - ints) > 0)
    # Get the position and the change. 
    pos = np.argmax(diff)
    change = list(seq['seq'])[pos]

    # Add the position and change as x and y
    loop_muts.loc[loop_muts['mutant']==g, 'x'] = pos + 1
    loop_muts.loc[loop_muts['mutant']==g, 'y'] = nt_idx[change]


# Compute the fractional width of the hpd
wt_vals = [pcut_stats[pcut_stats['mutant']=='WT12rss'], 
           ploop_stats[ploop_stats['mutant']=='WT12rss']]

for i, mut in enumerate([cut_muts, loop_muts]):
    avg_err = np.mean([mut['hpd_max'] - mut['median'],
                  mut['median'] - mut['hpd_min']])

    frac_err = (avg_err) / mut['median'].values

    mut['frac_err'] = frac_err * 7 

    # Add colors
    mut.sort_values('relative_prob', inplace=True)
   
cut_muts['color'] = sns.color_palette('PuOr', 
                n_colors=len(cut_muts)).as_hex()
loop_muts['color'] = sns.color_palette('PRGn_r', 
                n_colors=len(loop_muts)).as_hex()


# %%
# Find which mutants weren't queried
x, y = [], []
for i in range(1, 29):
    for j in range(4):
        if len(cut_muts[(cut_muts['x'] == i) & (cut_muts['y']==j)]) == 0:
            x.append(i)
            y.append(j)

# %%
#  # Visualization
ref_seq = vdj.io.endogenous_seqs()['reference']

cut_ax = bokeh.plotting.figure(width=700, height=200, x_axis_label='wild-type sequence',
                           y_axis_label='mutation', x_range=[-1, 29],
                           y_range=[-1, 4], title='cutting probability relative to wild type')
loop_ax = bokeh.plotting.figure(width=700, height=200, x_axis_label='wild-type sequence',
                            y_axis_label='mutation', x_range=[-1, 29],
                            y_range=[-1, 4], title='looping probability relative to wild type')
cut_ax.harea(0, 7, 4, color='red')
# Set the tick labels
for ax in [cut_ax, loop_ax]:
    ax.xaxis.ticker = np.arange(1, 29)
    ax.yaxis.ticker = [0, 1, 2, 3]
    ax.xaxis.major_label_overrides = {str(i):a for i, a in zip(np.arange(1, 29), list(ref_seq[0]))}
    ax.yaxis.major_label_overrides = {v:str(i) for i, v in nt_idx.items()} 
    ax.x(x, y, size=7, fill_color='lightgrey', line_color='black', 
          line_width=0.25, alpha=0.5)

    ax.square(np.arange(1, 29), ref_seq[1], size=10, fill_color='white', line_color='tomato')

cut_vals = cut_ax.circle(x='x', y='y', size='frac_err', line_width=0.5,  fill_color='color', 
        line_color='black', source=cut_muts)
loop_vals = loop_ax.circle(x='x', y='y', size='frac_err', line_width=0.5,  fill_color='color', 
        line_color='black', source=loop_muts)

# Plot the unexplored mutants
cut_hover = bokeh.models.HoverTool(renderers=[cut_vals], 
        tooltips=[('mutant', '@mutant'), ('mean', '@mean'), 
                  ('median', '@median'), ('mode','@mode'), 
                  ('minium 95% CR', '@hpd_min'), 
                  ('maximum 95% CR', '@hpd_max'),
                  ('change in probability', '@relative_prob')])
loop_hover = bokeh.models.HoverTool(renderers=[loop_vals], 
        tooltips=[('mutant', '@mutant'), ('mean', '@mean'), 
                  ('median', '@median'), ('mode','@mode'), 
                  ('minium 95% CR', '@hpd_min'), 
                  ('maximum 95% CR', '@hpd_max'),
                  ('change in probability', '@relative_prob')])

# highlight the regions
cut_ax.add_tools(cut_hover)
loop_ax.add_tools(loop_hover)
layout = bokeh.layouts.column(cut_ax, loop_ax)
# bokeh.io.show(layout)
# bokeh.io.output_file('./delta_p_cut.html')
# bokeh.io.save(layout)
bokeh.io.show(layout)


#%%


#%%

