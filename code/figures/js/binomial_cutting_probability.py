#-*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
from bokeh.models import ColumnDataSource, CustomJS,  HoverTool
from bokeh.models.glyphs import Patch
import bokeh.palettes
import bokeh.transform
import vdj.io
import vdj.viz
vdj.viz.plotting_style_bokeh()
bokeh.io.output_notebook()

#%%
# Load the summarized data and the posteriors
data = pd.read_csv('../../../data/pooled_cutting_probability.csv')
data = data[(data['hmgb1'] == 80) & (data['salt']=='Mg')]
posts = pd.read_csv('../../../data/pooled_cutting_probability_posteriors.csv')
posts = posts[(posts['hmgb1'] == 80) & (posts['salt']=='Mg')]


# Isoate only the point mutants for now
point = data[(data['n_muts'] == 1.0) & (data['mutant'] != 'V4-55')].copy()
point_posts = posts[(posts['n_muts'] == 1.0) & (posts['mutant'] != 'V4-55')].copy()

# Get the wild-type sequence information to generate the bubble plot
seqs = vdj.io.endogenous_seqs()
ref_seq = seqs['WT12rss'][0]
ref_idx = seqs['WT12rss'][1]

# Insert the x and y locations for the point mutants
for g, d in point.groupby('mutant'):
    seq = vdj.io.mutation_parser(g)
    seq_idx = seq['seq_idx']
    loc = np.argmax(seq_idx != ref_idx)
    val = seq_idx[loc]
    point.loc[point['mutant']==g, 'x'] = loc + 1
    point.loc[point['mutant']==g, 'y'] = val + 1

# Figure out the cutting probability of the wildtype
wt_prob = data[data['mutant']=='WT12rss']['mode'].values[0]
point['diff'] = point['mode'].values - wt_prob
point['size'] =  8 + (point['n_beads'] / 25)

# Define the wild-type data set
wt_data = data[data['mutant']=='WT12rss']
size = 8 + (wt_data['n_beads'].values[0] / 25)
_wt = pd.DataFrame(np.array([np.arange(len(ref_idx)) + 1, ref_idx + 1]).T, 
                    columns=['x', 'y'])
_wt['mode'] = wt_prob
_wt['size'] = size
_wt['diff'] = 0 
_wt['mutant'] = 'WT'
_wt['n_beads'] = wt_data['n_beads'].values[0]
_wt['n_cuts'] = wt_data['n_cuts'].values[0]
_wt['std'] = wt_data['std'].values[0]
# Determine where to plot the x's for unexplored mutations
xs, ys = [], []
for i in range(1, 29):
    for j in range(1, 5):
        if (len(point[(point['x']==i) & (point['y']==j)]) == 0) & (ref_idx[i - 1] != j - 1):
            xs.append(i)
            ys.append(j)


# Set the colors 
prob_cmap = bokeh.transform.linear_cmap(field_name='mode', palette=bokeh.palettes.Viridis[11], low=0, high=1)
diff_cmap = bokeh.transform.linear_cmap(field_name='diff', palette=bokeh.palettes.RdBu[11], low=-0.5, high=0.5)

# Define column data sources
muts = ColumnDataSource(point)
muts_post = ColumnDataSource(point_posts)
display = ColumnDataSource(dict(x=[], y=[], legend=[]))

# Write a callback to display the posterior of the hovered mutant
cb = bokeh.models.CustomJS(args=dict(mut_source=muts, post_source=muts_post, 
                                      post_source_display=display), code = """
        var data = post_source.data;
        var display = post_source_display.data;
        var mut_ind = cb_data.index['1d'].indices[0];
        var mut = mut_source.data['mutant'][mut_ind];
        var start = data['mutant'].indexOf(mut);
        var stop = data['mutant'].lastIndexOf(mut);
        console.log(mut)
        display['x'] = data['probability'].slice(start, stop+1);

        display['y'] = data['posterior'].slice(start, stop+1);
        display['legend'] = data['mutant'].slice(start, stop + 1);
        post_source_display.change.emit();
         """)


#%%
tooltips = [('Mutant', '@mutant'), ('# Beads', '@n_beads'), ('# Cuts', '@n_cuts'),
            ('cutting probability', '@mode'), ('standard deviation', '@std'),
            ('change from wild-type', '@diff')]

# Set up the figure canvas
prob_ax = bokeh.plotting.figure(width=600, height=200, 
                          x_axis_label='reference sequence',
                          y_axis_label='mutation',
                          x_range=[0, 29],
                          y_range=[0, 5],
                          title='cutting probability',
                          tools=[])
diff_ax = bokeh.plotting.figure(width=600, height=200, 
                          x_axis_label='reference sequence',
                          y_axis_label='mutation',
                          x_range=[0, 29],
                          y_range=[0, 5],
                          title='change in cutting probability', 
                          tools=[])
post_ax = bokeh.plotting.figure(width=600, height=300, x_axis_label='cutting probability',
                                y_axis_label='probability', 
                                title='posterior probability distribution',
                                y_range=[0, 0.08])



# Plot the mutants
prob_vals = prob_ax.circle(x='x', y='y', size='size', source=muts, fill_color=prob_cmap)
diff_vals = diff_ax.circle(x='x', y='y', size='size', source=muts, fill_color=diff_cmap)

# Plot the wild-type value
prob_ax.triangle(x='x', y='y', size='size', source=_wt, fill_color=prob_cmap,
              line_color='tomato', line_width=1.5, alpha=0.5)
diff_ax.triangle(x='x', y='y', size='size', source=_wt, fill_color=diff_cmap,
              line_color='dodgerblue', line_width=1.5, alpha=0.5)

# Plot the unexplored mutations
prob_ax.x(x=xs, y=ys, color='gray', alpha=0.5, size=8)
diff_ax.x(x=xs, y=ys, color='gray', alpha=0.5, size=8)

prob_hover = HoverTool(renderers= [prob_vals],
            tooltips=tooltips, callback=cb)
diff_hover = HoverTool(renderers= [diff_vals],
            tooltips=tooltips, callback=cb)
prob_hover.callback=cb
diff_hover.callback=cb
prob_ax.add_tools(prob_hover)
diff_ax.add_tools(diff_hover)

# Plot the wild-type posterior pdf
wt_post = posts[posts['mutant']=='WT12rss']
post_ax.line(x='probability', y='posterior', color='dodgerblue', legend='reference',
             source=wt_post)
post_ax.patch(x='probability', y='posterior', color='dodgerblue', source=wt_post, alpha=0.5, line_width=2)

# Plot the mutant posterior pdf
post_ax.line(x='x', y='y', color='tomato',  legend='legend', source=display, line_width=2)
post_ax.patch(x='x', y='y', color='tomato', source=display, alpha=0.5, line_width=2)
for p in [prob_ax, diff_ax]:
    p.yaxis.ticker = [1, 2, 3, 4]
    p.yaxis.major_label_overrides = {1:'A', 2:'C', 3:'G', 4:'T'}
    p.xaxis.ticker = np.arange(1, 29, 1)
    p.xaxis.major_label_overrides = {i+1:b for i, b in enumerate(list(ref_seq))}

col = bokeh.layouts.column(prob_ax, diff_ax, post_ax)
# bokeh.io.show(col)
bokeh.io.save(col, './cutting_prob.html')
#%%


#%wt_%
