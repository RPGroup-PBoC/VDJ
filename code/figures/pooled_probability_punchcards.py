# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.palettes
import vdj.io
import bokeh.transform
bokeh.io.output_notebook()

# Load teh cutting data
data = pd.read_csv('../../data/pooled_cutting_probability_summary.csv')
data.head()

#%% 
# Go through each mutant and compute the pooled probability
ref = vdj.io.endogenous_seqs()['reference']
ref_seq = ref[0]
ref = ref[1]
nt_idx = vdj.io.nucleotide_idx()
dfs = []
for g, d in data.groupby(['mutant']):
    # Compute the number of base pair differences and determine the x,y position
    seq = vdj.io.mutation_parser(g)
    n_diff = np.sum(np.abs(ref != seq['seq_idx']))
    if n_diff == 1:
        x = np.argmax(np.abs(ref - seq['seq_idx']))
        y = nt_idx[list(seq['seq'])[x]]
        x += 1
    else:
        x = None
        y = None
    df =  d.copy()
    df['mutant'] = g
    df['x'] = x
    df['y'] = y
    df['n_muts'] = n_diff
    dfs.append(df)
df = pd.concat(dfs)

#%%
# Compute the relative change in p_cut
ref_pcut = df[df['mutant']=='WT12rss']['median'].values[0]
df['relative_prob'] = df['median'] - ref_pcut
df['size'] = 5 / np.abs(df['hpd_min'] - df['hpd_max']) 
# Isolate the point mutants
point_muts = df[df['n_muts']==1]

# %%
# Find which mutants weren't queried
x, y = [], []
for i in range(1, 29):
    for j in range(4):
        if len(point_muts[(point_muts['x'] == i) & (point_muts['y']==j)]) == 0:
            x.append(i)
            y.append(j)
# %% 
# Define teh colors and assign to pcut values
palette = bokeh.palettes.PuOr11
colors = bokeh.transform.linear_cmap(field_name='relative_prob', low=-0.5, 
                                    high=0.5, low_color='purple', high_color='orange',
                                    palette=palette)
#%%
cut_ax = bokeh.plotting.figure(width=700, height=200, x_axis_label='wild-type sequence',
                           y_axis_label='mutation', x_range=[-1, 29],
                           y_range=[-1, 4], title='cutting probability relative to wild type')
# loop_ax = bokeh.plotting.figure(width=700, height=200, x_axis_label='wild-type sequence',
                            # y_axis_label='mutation', x_range=[-1, 29],
                            # y_range=[-1, 4], title='looping probability relative to wild type')
# cut_ax.harea(0, 7, 4, color='red')
# Set the tick labels

cut_ax.xaxis.ticker = np.arange(1, 29)
cut_ax.yaxis.ticker = [0, 1, 2, 3]
cut_ax.xaxis.major_label_overrides = {str(i):a for i, a in zip(np.arange(1, 29), list(ref_seq))}
cut_ax.yaxis.major_label_overrides = {'0':'A', '1':'C', '2':'G', '3':'T'} 
cut_ax.x(x, y, size=7, fill_color='lightgrey', line_color='black', 
          line_width=0.25, alpha=0.5)

cut_vals = cut_ax.circle(x='x', y='y', size='size', line_width=0.5,  fill_color=colors, 
        line_color='black', source=point_muts)
# loop_vals = loop_ax.circle(x='x', y='y', size='frac_err', line_width=0.5,  fill_color='color', 
#         line_color='black', source=loop_muts)

# Plot the unexplored mutants
cut_hover = bokeh.models.HoverTool(renderers=[cut_vals], 
        tooltips=[('mutant', '@mutant'), ('median', '@median'), 
                  ('mean', '@mean'), ('mode', '@mode'), 
                  ('change in probability', '@relative_prob')])
#loop_hover = bokeh.models.HoverTool(renderers=[loop_vals], 
#         tooltips=[('mutant', '@mutant'), ('mean', '@mean'), 
#                   ('median', '@median'), ('mode','@mode'), 
#                   ('minium 95% CR', '@hpd_min'), 
#                   ('maximum 95% CR', '@hpd_max'),
#                   ('change in probability', '@relative_prob')])

# highlight the regions
cut_ax.add_tools(cut_hover)
# loop_ax.add_tools(loop_hover)
# layout = bokeh.layouts.column([cut_ax])
# bokeh.io.show(layout)
# bokeh.io.output_file('./delta_p_cut.html')
# bokeh.io.save(layout)
bokeh.io.show(cut_ax)


#%%



#%%
