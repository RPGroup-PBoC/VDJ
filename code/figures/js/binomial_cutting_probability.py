#-*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
from bokeh.models import ColumnDataSource, CustomJS,  HoverTool
from bokeh.models.glyphs import Patch
from bokeh.resources import INLINE
import bokeh.palettes
import bokeh.transform
import vdj.io
import vdj.viz
vdj.viz.plotting_style_bokeh()
bokeh.io.output_notebook(resources=INLINE)

#%%
# Load the summarized data and the posteriors
data = pd.read_csv('../../../data/pooled_cutting_probability.csv')
data = data[(data['hmgb1'] == 80) & (data['salt']=='Mg')]
posts = pd.read_csv('../../../data/pooled_cutting_probability_posteriors.csv')
posts = posts[(posts['hmgb1'] == 80) & (posts['salt']=='Mg')]
dwell = pd.read_csv('../../../data/compiled_dwell_times.csv')
dwell = dwell[(dwell['salt']=='Mg') & (dwell['hmgb1']==80)]

# Isoate only the point mutants for now
point = data[data['n_muts'] == 1.0].copy()
point_posts = posts[posts['n_muts'] == 1.0].copy()

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
    dwell.loc[dwell['mutant']==g, 'n_muts'] = seq['n_muts']

dwell = dwell[dwell['n_muts']==1]
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

# pre-bin the histogram
bins = np.linspace(0, dwell['dwell_time_min'].max(), 30)
hist_dfs = []
cut_dfs = []
for g, d in dwell.groupby(['mutant']):
    # Dist all dwell times
    hist, bins = np.histogram(d['dwell_time_min'], bins=bins)
    # hist = hist / np.sum(hist)
    _df = pd.DataFrame(np.array([hist, bins[:-1], bins[1:]]).T, columns=['top', 'left', 'right'])
    _df['mutant'] = g
    _df['bottom'] = 0
    hist_dfs.append(_df)

    # Hist only the cutting events
    _d = d[d['cut']==1]
    cut_hist, cut_bins = np.histogram(_d['dwell_time_min'], bins=bins)
    # cut_hist = cut_hist / np.sum(hist)
    _df = pd.DataFrame(np.array([cut_hist, cut_bins[:-1], cut_bins[1:]]).T, 
                       columns=['top', 'left', 'right'])
    _df['mutant'] = g
    _df['bottom'] = 0 
    cut_dfs.append(_df)
hist_df = pd.concat(hist_dfs)
cut_hist_df = pd.concat(cut_dfs)

# Generate the column data source
hist_source = ColumnDataSource(hist_df)
cut_source = ColumnDataSource(cut_hist_df)
wt_hist = hist_df[hist_df['mutant']=='WT12rss'] 
wt_cut_hist = cut_hist_df[cut_hist_df['mutant']=='WT12rss'] 
hist_display = ColumnDataSource(dict(left=[], right=[], bottom=[], top=[]))
cut_hist_display = ColumnDataSource(dict(left=[], right=[], bottom=[], top=[]))

# Define column data sources
muts = ColumnDataSource(point)
muts_post = ColumnDataSource(point_posts)
display = ColumnDataSource(dict(x=[], y=[], legend=[]))

# Write a callback to display the posterior of the hovered mutant
cb = bokeh.models.CustomJS(args=dict(mut_source=muts, post_source=muts_post, 
                                      post_source_display=display,
                                      hist_source=hist_source, 
                                      hist_source_display=hist_display,
                                      cut_source=cut_source,
                                      cut_hist_display=cut_hist_display), code = """
        var data = post_source.data;
        var hist_data = hist_source.data;
        var cut_data = cut_source.data;
        var display = post_source_display.data;
        var hist_display = hist_source_display.data;
        var cut_display = cut_hist_display.data;
        var mut_ind = cb_data.index['1d'].indices[0];
        var mut = mut_source.data['mutant'][mut_ind];
        var start = data['mutant'].indexOf(mut);
        var stop = data['mutant'].lastIndexOf(mut);
        var hist_start = hist_data['mutant'].indexOf(mut);
        var hist_stop = hist_data['mutant'].lastIndexOf(mut);
        var cut_start = cut_data['mutant'].indexOf(mut);
        var cut_stop = cut_data['mutant'].lastIndexOf(mut);
        display['x'] = data['probability'].slice(start, stop+1);
        display['y'] = data['posterior'].slice(start, stop+1);
        display['legend'] = data['mutant'].slice(start, stop + 1);
        hist_display['legend'] = hist_data['mutant'].slice(hist_start, hist_stop+1)
        hist_display['left'] = hist_data['left'].slice(hist_start, hist_stop+1);
        hist_display['right'] = hist_data['right'].slice(hist_start, hist_stop+1);
        hist_display['top'] = hist_data['top'].slice(hist_start, hist_stop+1);
        hist_display['bottom'] = hist_data['bottom'].slice(hist_start, hist_stop+1);
        cut_display['left'] = cut_data['left'].slice(cut_start, cut_stop+1);
        cut_display['right'] = cut_data['right'].slice(cut_start, cut_stop+1);
        cut_display['top'] = cut_data['top'].slice(cut_start, cut_stop+1);
        cut_display['bottom'] = cut_data['bottom'].slice(cut_start, cut_stop+1);
        cut_display['legend'] = cut_data['mutant'].slice(cut_start, cut_stop+1)
        post_source_display.change.emit();
        cut_hist_display.change.emit();
        hist_source_display.change.emit();
         """)


#%%
tooltips = [('Mutant', '@mutant'), ('# Loops', '@n_beads'), ('# Cuts', '@n_cuts'),
            ('cutting probability', '@mode'), ('standard deviation', '@std'),
            ('change from reference', '@diff')]

# Set up the figure canvas
diff_ax = bokeh.plotting.figure(width=700, height=200, 
                          x_axis_label='reference sequence',
                          y_axis_label='mutation',
                          x_range=[0, 29],
                          y_range=[0, 5],
                          title='change in cutting probability', 
                          tools=[])
post_ax = bokeh.plotting.figure(width=350, height=300, x_axis_label='cutting probability',
                                y_axis_label='probability', 
                                title='posterior probability distribution',
                                y_range=[0, 0.08])
dist_ax = bokeh.plotting.figure(width=350, height=300,
                          x_axis_label='paired-complex dwell time [min]',
                                y_axis_label='number of loops', 
                                title='dwell time distributions')



# dist_ax.quad(bottom='bottom', top='top', left='left', right='right',
#     source=wt_hist, fill_alpha=0.2, line_width=1, legend='reference',
#     line_color='dodgerblue', fill_color='dodgerblue')
# dist_ax.quad(bottom='bottom', top='top', left='left', right='right',
#     source=wt_cut_hist, fill_alpha=0.4, line_width=1, legend='reference (cut)',
#     line_color='grey', fill_color='grey')

dist_ax.quad(bottom='bottom', top='top', left='left', right='right',
    source=hist_display, fill_alpha=0.4, line_width=1, legend='all loops',
    line_color='tomato', fill_color='tomato')
dist_ax.quad(bottom='bottom', top='top', left='left', right='right',
    source=cut_hist_display, fill_alpha=0.0, line_width=1, legend='cut loops',
    line_color='black', fill_color='black', hatch_pattern='/', hatch_color='black')


# Plot the mutants
diff_vals = diff_ax.circle(x='x', y='y', size='size', source=muts, fill_color=diff_cmap)

# Plot the wild-type value
diff_ax.triangle(x='x', y='y', size='size', source=_wt, fill_color=diff_cmap,
              line_color='dodgerblue', line_width=1.5, alpha=0.5)

# Plot the unexplored mutations
diff_ax.x(x=xs, y=ys, color='gray', alpha=0.5, size=8)

diff_hover = HoverTool(renderers= [diff_vals],
            tooltips=tooltips, callback=cb)

diff_hover.callback=cb
diff_ax.add_tools(diff_hover)

# Plot the wild-type posterior pdf
wt_post = posts[posts['mutant']=='WT12rss']
post_ax.line(x='probability', y='posterior', color='dodgerblue', legend='reference',
             source=wt_post)
post_ax.patch(x='probability', y='posterior', color='dodgerblue', source=wt_post, alpha=0.5, line_width=2)

# Plot the mutant posterior pdf
post_ax.line(x='x', y='y', color='tomato',  legend='legend', source=display, line_width=2)
post_ax.patch(x='x', y='y', color='tomato', source=display, alpha=0.5, line_width=2)
diff_ax.yaxis.ticker = [1, 2, 3, 4]
diff_ax.yaxis.major_label_overrides = {1:'A', 2:'C', 3:'G', 4:'T'}
diff_ax.xaxis.ticker = np.arange(1, 29, 1)
diff_ax.xaxis.major_label_overrides = {i+1:b for i, b in enumerate(list(ref_seq))}

row = bokeh.layouts.row(post_ax, dist_ax)
col = bokeh.layouts.column(diff_ax, row)
bokeh.io.show(col)
bokeh.io.save(col, './cutting_prob.html')
#%%


#%wt_%
