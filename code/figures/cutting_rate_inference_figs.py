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
samples = pd.read_csv('../../data/pooled_cutting_rate_samples.csv')
stats = pd.read_csv('../../data/pooled_cutting_rate.csv')

# Consider only the  single mutants
points = stats[(stats['n_muts']==1) & (stats['parameter']=='tau')].copy()

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
    points.loc[points['mutant']==m, 'size'] = 0.1 * (_m['hpd_max'] - _m['hpd_min'].values) 

# Bin the sample data
binned_dfs = []
bins = np.linspace(0, 1110, 100)
for g, d in samples.groupby('mutant'):
    hist, bins = np.histogram(d['tau'], bins)
    df = pd.DataFrame([])
    df['bins'] = bins[:-1]
    df['hist'] = hist / np.sum(hist)
    df['median'] = np.median(d['tau'])
    df['mutant'] =  g
    df['alpha'] = 0
    binned_dfs.append(df)
binned = pd.concat(binned_dfs)

# # Define the data sets as sources
dwell_source = bokeh.models.ColumnDataSource(data)
stats_source = bokeh.models.ColumnDataSource(points)
bin_source = bokeh.models.ColumnDataSource(binned)
bin_display = bokeh.models.ColumnDataSource({'bins': [0], 'hist':[0]})

# # Define the dropdown menu
# drop = bokeh.models.widgets.Select(title='Select Mutant', value='WT12rss', 
#                                     options=list(stats['mutant'].unique()))


# %%
greens = bokeh.palettes.Blues9
tau = stats[stats['parameter']=='tau']['median']
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
dist_ax = bokeh.plotting.figure(width=canvas_width, height=200, x_range=[1, 1150],
                               x_axis_label='τ [sec/cut]',
                               y_axis_label='probability')

# Plot the  punch cards
rate_vals = rate_ax.circle(x='x', y='y', fill_color=rate_colors, source=stats_source, 
                line_color='black', size='size', hover_fill_color='tomato')


# Plot the unexplored mutations
x, y = [], [] 
for _x in range(1, 29):
    for _y in range(1, 5):
        if len(points[(points['x']==_x) & (points['y']==_y)]) == 0:
            x.append(_x)
            y.append(_y)
rate_ax.x(x=x, y=y, color='grey', alpha=0.5, size=8)

# Plot the wild-type values
_wt = stats[(stats['parameter']=='tau') & (stats['mutant']=='WT12rss')]
wt = pd.DataFrame([])
wt['x'] = np.arange(1, 29)
wt['y'] = ref_idx + 1
wt['mutant'] = 'WT'
wt['mean'] = _wt['mean'].values[0]
wt['mode'] = _wt['mode'].values[0]
wt['hpd_max'] = _wt['hpd_max'].values[0]
wt['hpd_min'] = _wt['hpd_min'].values[0]
wt['median'] = _wt['median'].values[0]
wt['size'] = 6 * (wt['median'] / (_wt['hpd_max'].values[0] - _wt['hpd_min'].values[0]))
wt_rate_vals = rate_ax.circle(x='x', y='y', fill_color=rate_colors, source=wt, 
            line_color='tomato') 


# Bin and show the distribution for the wild-type
wt_samps = binned[binned['mutant']=='WT12rss']
dist_ax.step(wt_samps['bins'], wt_samps['hist'], line_width=1, color='blue', legend='wild type')
dist_ax.step(x='bins', y='hist', line_width=1, color='tomato', legend='legend', 
            source=bin_display)


# Add hover tools
cb = bokeh.models.CustomJS(args=dict(mut_source=stats_source, bin_source=bin_source, 
                                      bin_source_display=bin_display), code = """
        var data = bin_source.data;
        var display = bin_source_display.data;
        var mut_ind = cb_data.index['1d'].indices[0];
        var mut = mut_source.data['mutant'][mut_ind];
        var start = data['mutant'].indexOf(mut);
        var stop = data['mutant'].lastIndexOf(mut);
        display['bins'] = data['bins'].slice(start, stop+1);
        display['hist'] = data['hist'].slice(start, stop+1);
        console.log(display)
        bin_source_display.change.emit();
         """)

rate_hover = bokeh.models.HoverTool(renderers=[rate_vals], 
        tooltips=[('mutant', '@mutant'), ('median τ [sec/cut]', '@median'), 
                  ('mean τ [sec/cut]', '@mean'), ('most probable τ [sec/cut] ', '@mode'), 
                  ('maximum 95% CR [sec/cut]', '@hpd_max'), 
                  ('minimum 95% CR [sec/cut]', '@hpd_min'),
                  ('index', '$index')], callback=cb)


dist_ax.step(x='bins', y='hist', color='tomato', line_width=2, source=bin_display)
rate_hover.callback = cb
wt_rate_hover = bokeh.models.HoverTool(renderers=[wt_rate_vals], 
        tooltips=[('mutant', '@mutant'), ('median τ [sec/cut]', '@median'), 
                  ('mean τ [sec/cut]', '@mean'), ('most probable τ [sec/cut] ', '@mode'), 
                  ('maximum 95% CR [sec/cut]', '@hpd_max'), 
                  ('minimum 95% CR [sec/cut]', '@hpd_min')])
rate_ax.add_tools(rate_hover)
rate_ax.add_tools(wt_rate_hover)




# Define the layout

col = bokeh.layouts.column(dwell_ax, rate_ax, rel_rate_ax, dist_ax)

bokeh.io.show(col)


bokeh.io.save(col, './test.html')

#%%
p = bokeh.plotting.figure(width=900)
bins = np.linspace(0, 1000, 50)
hist, bins = np.histogram(samples[samples['mutant']=='WT12rss']['tau'], bins=bins)
p.step(bins[:-1], hist, legend='WT', color='black')
hist, bins = np.histogram(samples[samples['mutant']=='12SpacG11T']['tau'], bins=bins)
p.step(bins[:-1], hist, legend='12SpacG11T', color='tomato')
hist, bins = np.histogram(samples[samples['mutant']=='12HeptT6A']['tau'], bins=bins)
p.step(bins[:-1], hist, legend='12HeptG6A', color='dodgerblue')
hist, bins = np.histogram(samples[samples['mutant']=='12SpacC4T']['tau'], bins=bins)
p.step(bins[:-1], hist, legend='12SpacC4T', color='purple')
hist, bins = np.histogram(samples[samples['mutant']=='12HeptC3G']['tau'], bins=bins)
p.step(bins[:-1], hist, legend='12HeptC3G', color='orange')


bokeh.io.show(p)

#%%
samples[samples['mutant']=='12SpacG4T']

#%%
