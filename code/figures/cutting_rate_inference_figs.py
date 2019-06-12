# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.models
import bokeh.layouts
import bokeh.colors
import bokeh.transform
import scipy.stats
import vdj.io
bokeh.io.output_notebook()

# Load the statistics
data = pd.read_csv('../../data/compiled_dwell_times.csv')
posteriors = pd.read_csv('../../data/cutting_rate_posteriors.csv')
stats = pd.read_csv('../../data/cutting_rate_analytic_summary.csv')

# Subsample the posteriors
# posteriors = posteriors[posteriors['tau'] < 1200]

# %%
# Consider only the single mutants
points = stats[(stats['n_muts']<=1)].copy()

# Load the reference sequence and define the x, y positions of the base changes
seqs = vdj.io.endogenous_seqs()
nt_idx = vdj.io.nucleotide_idx()
ref = seqs['reference']
ref_idx = ref[1] 
for m in points['mutant'].unique():
    mut = vdj.io.mutation_parser(m)
    mut_idx = mut['seq_idx']

    # Find the location of the disagreement
    ind = np.argmax(ref_idx != mut_idx) 
    x = ind + 1
    y = nt_idx[mut['seq'][ind]]
    _m = points[points['mutant']==m]
    # Add it to the data frame
    points.loc[points['mutant']==m, 'x'] = x
    points.loc[points['mutant']==m, 'y'] = y + 1
    points.loc[points['mutant']==m, 'base'] = mut['seq'][ind]
    points.loc[points['mutant']==m, 'size'] = np.log10(100/_m['width'])


# Compute the CDFS for each
cdf_dfs = []
for g, d in data.groupby(['mutant']):
        _d = d[d['cut']==1]
        x = np.sort(_d['dwell_time_min'] - 0.35)
        y = np.arange(1, len(_d) + 1) / len(_d)
        _df = pd.DataFrame(np.array([x, y]).T, columns=['x', 'y'])
        _df['mutant'] = g
        if g in points['mutant'].unique():
                cdf_dfs.append(_df)
cdf_df = pd.concat(cdf_dfs)


# %%
# # Define the data sets as sources
dwell_source = bokeh.models.ColumnDataSource(data)
stats_source = bokeh.models.ColumnDataSource(points)
post_source = bokeh.models.ColumnDataSource(posteriors[posteriors['n_muts']==1])
post_display = bokeh.models.ColumnDataSource({'x':[], 
                                                'y':[], 'legend':[]})

# # Define the dropdown menu
# drop = bokeh.models.widgets.Select(title='Select Mutant', value='WT12rss', 
#                                     options=list(stats['mutant'].unique()))


# %%
greens = bokeh.palettes.Blues9
tau = stats['median']
rate_colors = bokeh.transform.linear_cmap(palette=greens, 
        low = tau.min(), high=tau.max(), 
        low_color=greens[0], high_color=greens[-1],
        field_name='median')

# %%
# Set up the figures
canvas_width = 800
dwell_ax = bokeh.plotting.figure(width=canvas_width, height=400, 
                                x_axis_label='dwell time [s]',
                                y_axis_label='cumulative distribution')
rate_ax = bokeh.plotting.figure(width=canvas_width, height=200, 
                               x_range=[0, 29], y_range=[0,5],
                               x_axis_label='reference sequence', 
                               y_axis_label='mutation')
rel_rate_ax = bokeh.plotting.figure(width=canvas_width, height=200,
                               x_range=[1, 29], y_range=[1, 4],
                               x_axis_label='reference sequence',
                               y_axis_label='mutation')
dist_ax = bokeh.plotting.figure(width=canvas_width, height=200, x_range=[1, 60],
                               x_axis_label='τ [min/cut]',
                               y_axis_label='probability')

# Plot the  punch cards
rate_vals = rate_ax.circle(x='x', y='y', fill_color=rate_colors, source=stats_source, 
                line_color='black', size=20, hover_fill_color='tomato')


# Plot the unexplored mutations
x, y = [], [] 
for _x in range(1, 29):
    for _y in range(1, 5):
        if len(points[(points['x']==_x) & (points['y']==_y)]) == 0:
            x.append(_x)
            y.append(_y)
rate_ax.x(x=x, y=y, color='grey', alpha=0.5, size=8)

# Plot the wild-type values
_wt = stats[stats['mutant']=='12SpacG11T']
wt = pd.DataFrame([])
wt['x'] = np.arange(1, 29)
wt['y'] = ref_idx + 1
wt['mutant'] = 'WT'
wt['95_low'] = _wt['95_low'].values[0] 
wt['95_high'] = _wt['95_high'].values[0]
wt['median'] = _wt['median'].values[0]
wt['size'] = 30 * np.log10(500 / _wt['width'])
wt_rate_vals = rate_ax.circle(x='x', y='y', fill_color=rate_colors, source=wt, 
            line_color='tomato', size=20) 


#  Pot the wild-type dwell times
_wt_dwell = cdf_df[cdf_df['mutant']=='12SpacG11T']
dwell_ax.step(x='x', y='y', line_width=1, color='blue', source=_wt_dwell)
dwell_ax.circle(x='x', y='y', line_width=1, line_color='blue', source=_wt_dwell,
                fill_color='white')

# Plot the theoretical cdf
time = np.linspace(0, 50, 100)
wt_tau = stats[stats['mutant']=='12SpacG11T']['median']
theo_cdf = 1 - np.exp(-(time - 0.35)/wt_tau.values[0])
dwell_ax.step(time, theo_cdf, line_width=3, color='blue', alpha=0.5)
cut = data[(data['mutant']=='12SpacG11T') & (data['cut']==1)]['dwell_time_min']






bokeh.io.show(dwell_ax)

# %%
# Bin and show the distribution for the wild-type
wt_samps = posteriors[posteriors['mutant']=='WT12rss']
dist_ax.line(x='tau', y='posterior_pdf', line_width=1, color='blue', legend='wild type', source=wt_samps)

dist_ax.line(x='x', y='y', color='tomato', line_width=2, source=post_display)

# %%
# Add hover tool
cb = bokeh.models.CustomJS(args=dict(mut_source=stats_source, post_source=post_source, 
                                      post_source_display=post_display), code = """
        var data = post_source.data;
        var display = post_source_display.data;
        var mut_ind = cb_data.index['1d'].indices[0];
        var mut = mut_source.data['mutant'][mut_ind];
        var start = data['mutant'].indexOf(mut);
        var stop = data['mutant'].lastIndexOf(mut);
        display['x'] = data['tau'].slice(start, stop+1);
        display['y'] = data['posterior_pdf'].slice(start, stop+1);
        display['legend'] = data['mutant'].slice(start, stop + 1);
        post_source_display.change.emit();
         """)

rate_hover = bokeh.models.HoverTool(renderers=[rate_vals], 
        tooltips=[('mutant', '@mutant'), ('median τ [min/cut]', '@median'), 
                  ('maximum 95% CR [min/cut]', '@95_high'), 
                  ('minimum 95% CR [min/cut]', '@95_low'),
                  ('# loops', '@n_loops')], callback=cb)


rate_hover.callback = cb
wt_rate_hover = bokeh.models.HoverTool(renderers=[wt_rate_vals], 
        tooltips=[('mutant', '@mutant'), ('median τ [min/cut]', '@median'), 
                  ('maximum 95% CR [min/cut]', '@95_high'), 
                  ('minimum 95% CR [min/cut]', '@95_low'),
                  ('# loops', '@n_loops')])
rate_ax.add_tools(rate_hover)
rate_ax.add_tools(wt_rate_hover)

# Define the layout
col = bokeh.layouts.column(rate_ax, dist_ax)
bokeh.io.show(col)
bokeh.io.save(col, './cutting_rate.html')

#%%


#%%
