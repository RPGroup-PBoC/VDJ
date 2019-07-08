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

# Rename some of the mutant names to correct endogenous names
endog_names = {'WT12rss' : 'V4-57-1 (ref)',
                'DFL161' : 'DFL 16.1-5\'',
                'DFL1613' : 'DFL 16.1-3\'',
                '12SpacC1A' : 'V4-55'}

# Isolate only the point mutants for now
endog = data[(~data['mutant'].str.startswith('12')) | (data['mutant']=='12SpacC1A')].copy()
endog_posts = posts[(~posts['mutant'].str.startswith('12')) | (data['mutant']=='12SpacC1A')].copy()

endog = endog.replace({'mutant' : endog_names})
endog_posts = endog_posts.replace({'mutant' : endog_names})

endog_ordering_y = {'DFL 16.1-3\'' : 11.0,
                    'DFL 16.1-5\'' : 10.0,
                    'V1-135' : 9.0,
                    'V9-120' : 8.0,
                    'V10-96' : 7.0,
                    'V19-93' : 6.0,
                    'V4-57-1 (ref)' : 5.0,
                    'V4-55' : 4.0,
                    'V5-43' : 3.0,
                    'V8-18' : 2.0,
                    'V6-17' : 1.0,
                    'V6-15' : 0.0
                    }
#%%

# Insert the y locations for the point mutants
for g, d in endog.groupby('mutant'):
    val = endog_ordering_y[g]
    endog.loc[endog['mutant']==g, 'y'] = val + 1

#%%
dwell = dwell[(~dwell['mutant'].str.startswith('12')) | (dwell['mutant']=='12SpacC1A')]
dwell = dwell.replace({'mutant' : endog_names})
data = data.replace({'mutant' : endog_names})
# Figure out the cutting probability of the wildtype
wt_prob = data[data['mutant']=='V4-57-1 (ref)']['mode'].values[0]
endog['diff'] = endog['mode'].values - wt_prob
endog['size'] =  8 + (endog['n_beads'] / 25)

# Define the wild-type data set
wt_data = data[data['mutant']=='V4-57-1 (ref)']
size = 8 + (wt_data['n_beads'].values[0] / 25)
_wt = pd.DataFrame(np.array([np.arange(len(ref_idx)) + 1, ref_idx + 1]).T, 
                    columns=['y'])
_wt['mode'] = wt_prob
_wt['size'] = size
_wt['diff'] = 0 
_wt['mutant'] = 'V4-57-1 (ref)'
_wt['n_beads'] = wt_data['n_beads'].values[0]
_wt['n_cuts'] = wt_data['n_cuts'].values[0]
_wt['std'] = wt_data['std'].values[0]


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
muts = ColumnDataSource(endog)
muts_post = ColumnDataSource(endog_posts)
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
                          y_range=[0, len(endog_ordering_y)],
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
wt_post = endog_posts[endog_posts['mutant']=='WT12rss']
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
