# -*- coding: utf-8 -*-
"""
Builds a Bokeh appelet for exploring the data and visualizing the inference
statistics.
"""
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, GroupFilter 
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.models.glyphs import HBar
from bokeh.embed import components
bokeh.plotting.output_file('./vdj_voyager.html')
import imp
import vdj.io
imp.reload(vdj.io)

# Load the datasets
dwell_times = pd.read_csv('../../../data/compiled_dwell_times.csv')
dwell_times = dwell_times[dwell_times['mutant'] != 'analysis']
f_looped = pd.read_csv('../../../data/compiled_looping_fraction.csv')
f_looped = f_looped[f_looped['mutant'] != 'analysis']
fates = pd.read_csv('../../../data/compiled_cutting_events.csv')
fates = fates[fates['mutant'] != 'analysis']

# Instantiate the dropdown menu
selector = Select(title='Mutant', value='WT12rss', 
                  options=list(dwell_times.mutant.unique()))

#Get the sequences and define the colors for all mutants
ref = vdj.io.endogenous_seqs()['reference'][0]
seq_dfs = []
for m in fates.mutant.unique():
    seq = vdj.io.mutation_parser(m)['seq']
    colors = []
    for r, s in zip(list(ref), list(seq)):
        if r == s:
            colors.append('black')
        else:
            colors.append('dodgerblue')
    x = np.arange(len(seq))
    y = np.zeros_like(x)
    _df = pd.DataFrame([])
    _df['x'] = x
    _df['y'] = y
    _df['sequence'] = list(seq)
    _df['colors'] = list(colors)
    _df['mutant'] = m
    seq_dfs.append(_df)
seq_df = pd.concat(seq_dfs)

# Compute the ECDFs for the dwell times
pooled_dwell_times= []
for g, d in dwell_times.groupby('mutant'):
    # Compute the pooled ecdf
    x = np.sort(d['dwell_time_ms'])
    y = np.linspace(0, 1, len(d))
    df = pd.DataFrame([])
    df['pooled_dwell_time'] = x
    df['ecdf'] = y
    df['mutant'] = g
    pooled_dwell_times.append(df)

# Compute the replicate ecdfs
rep_dwell_times = []
for g, d in dwell_times.groupby(['mutant', 'replicate']):
    x  = np.sort(d['dwell_time_ms'])
    y = np.linspace(0, 1, len(d))
    df = pd.DataFrame([])
    df['dwell'] = x
    df['ecdf'] = y
    df['mutant'] = g[0] 
    df['replicate'] = g[1]
    rep_dwell_times.append(df)

# Compute the total bead counts for that particular mutant
pooled_fates = pd.DataFrame([])
for g, d in fates.groupby(['mutant']):
    pooled_fates = pooled_fates.append({'fate':'cut','value':d['n_cuts'].values.sum(),
                        'mutant':g}, ignore_index=True)
    pooled_fates = pooled_fates.append({'fate':'uncut','value':d['n_beads'].values.sum() -\
                        d['n_cuts'].values.sum(), 'mutant':g},
                        ignore_index=True)

# Compute the pooled probability of cutting
cut_prob = pd.DataFrame()
for g, d in fates.groupby(['mutant']):
    p_cut = d['n_cuts'].sum() / d['n_beads'].sum()
    cut_err = (d['n_cuts'].sum() * (d['n_beads'].sum() - d['n_cuts'].sum())) /\
                (d['n_beads'].sum()**3)
    low = p_cut - cut_err
    high = p_cut + cut_err
    val = {'mutant':g, 'p_cut':p_cut, 
            'p_cut_err':cut_err, 'low':low, 'high':high, 'y': 0}
    cut_prob = cut_prob.append(val, ignore_index=True)


# Compute the replicate level probability of cutting
rep_cut_prob = pd.DataFrame()
for g, d in fates.groupby(['mutant', 'replicate']):
    p_cut = d['n_cuts'].values[0] / d['n_beads'].values[0]
    cut_err = (d['n_cuts'].values[0] * (d['n_beads'].values[0] - d['n_cuts'].values[0])) / (d['n_beads'].values[0]**3)
    low = p_cut - cut_err
    high = p_cut + cut_err
    val = {'mutant':g[0], 'p_cut':p_cut, 
            'p_cut_err':cut_err, 'low':low, 'high':high, 'y': np.random.normal(0, 0.01)}
    rep_cut_prob = rep_cut_prob.append(val, ignore_index=True)
rep_cut_prob.dropna(inplace=True)


# Compute the fractional time in looped state for all data.
floop_pooled = f_looped.groupby(['mutant'])['fraction_looped'].agg(('mean', 'sem')).reset_index()
floop_pooled['y'] = 0
floop_rep = f_looped.groupby(['mutant', 'replicate'])['fraction_looped'].agg(('mean', 'sem')).reset_index()
floop_rep['y'] = np.random.normal(0, 0.01, len(floop_rep))

# Concatenate the dwell time distribution dataframes and generate CDS
rep_dwell = pd.concat(rep_dwell_times)
pooled_dwell = pd.concat(pooled_dwell_times)
rep_source = ColumnDataSource(rep_dwell)
pooled_source = ColumnDataSource(pooled_dwell)
fate_source = ColumnDataSource(pooled_fates)
pcut_source = ColumnDataSource(cut_prob)
rep_pcut_source = ColumnDataSource(rep_cut_prob)
floop_pooled_source = ColumnDataSource(floop_pooled)
floop_rep_source = ColumnDataSource(floop_rep)
seq_source = ColumnDataSource(seq_df)

# Define the filter for mutant selection
filter = GroupFilter(column_name="mutant", group="WT12rss")

# Assign the views
pooled_view = CDSView(source=pooled_source, filters=[filter])
rep_view = CDSView(source=rep_source, filters=[filter])
fate_view = CDSView(source=fate_source, filters=[filter] )
pcut_view = CDSView(source=pcut_source, filters=[filter])
rep_pcut_view = CDSView(source=rep_pcut_source, filters=[filter])
floop_view = CDSView(source=floop_pooled_source, filters=[filter])
floop_rep_view = CDSView(source=floop_rep_source, filters=[filter])
seq_view = CDSView(source=seq_source, filters=[filter])


# Define the callback args and callback code
cb_args = {'pv':pooled_view, 'rv':rep_view, 'fv':fate_view, 'pcv':pcut_view, 
           'rpcv':rep_pcut_view, 'flv':floop_view, 'flrv':floop_rep_view,
           'seqv':seq_view,
          'sel':selector, 'filter':filter,
          'ps':pooled_source, 'rs':rep_source, 'fs':fate_source, 
          'pcs':pcut_source, 'rpcs':rep_pcut_source, 'fls':floop_pooled_source,
          'flrs':floop_rep_source, 'seqs':seq_source}
dwell_cb = CustomJS(args=cb_args, 
    code="""
    var mut = sel.value
    filter.group = mut;
    pv.filters[0] = filter;
    seqv.filters[0] = filter;
    rv.filters[0] = filter;
    fv.filters[0] = filter;
    pcv.filters[0] = filter;
    rpcv.filters[0] = filter;
    flv.filters[0] = filter;
    flrv.filters[0] = filter;
    fls.data.view = pv;
    seqs.data.view = pv;
    flrs.data.view = pv;
    ps.data.view = pv;
    pcs.data.view = pv;
    rpcs.data.view = pv;
    rs.data.view = rv;
    fs.data.view = fv;
    fls.data.view = fv;
    flrs.data.view = fv;
    rs.change.emit();
    ps.change.emit();
    pcs.change.emit();
    rpcs.change.emit();
    fs.change.emit();
    fls.change.emit();
    flrs.change.emit();
    seqs.change.emit();
""")

# Dwell time figure canvas
dwell_ax = bokeh.plotting.figure(width=800, height=400, 
                           x_axis_label='dwell time [ms]',
                           y_axis_label='cumulative distribution',
                           title=f'paired complex dwell time',
                           x_range=[np.min(dwell_times['dwell_time_ms']), 
                                    np.max(dwell_times['dwell_time_ms'])])
_fates = ['uncut', 'cut']

# Cutting fraction canvas
cut_ax = bokeh.plotting.figure(width=200, height=150, 
                               x_axis_label='number of events',
                               y_range=_fates, title='bead fate')
                               
# Cutting Probability Canvas
pcut_ax = bokeh.plotting.figure(width=300, height=150, x_range=[-0.1, 1.1],
                               x_axis_label='probability', y_range=[-0.1, 0.1],
                title=f'empirical cutting probability')
pcut_ax.yaxis.visible = False

# Looping Fraction Canvas
floop_ax = bokeh.plotting.figure(width=300, height=150,
        title=f'fraction of time looped',
        x_axis_label='fractional time', y_range=[-0.1, 0.1], x_range=[-0.1, 0.5])
floop_ax.yaxis.visible = False

# Sequence canvas
seq = bokeh.plotting.figure(width=600, height=40, y_range=[0, 0.1], x_range=[0, 28])
seq.xaxis.visible=False
seq.yaxis.visible=False
glyph = bokeh.models.glyphs.Text(x='x', y='y', text='sequence', text_color='colors', 
        text_font_size='28px', text_font='Courier')
seq.add_glyph(seq_source, glyph, view=seq_view)

# Plot the dwell time distributions
dwell_ax.circle(x=[], y=[], line_width=2, color='slategrey', legend='replicate')
dwell_ax.circle(x='dwell', y='ecdf', source=rep_source, view=rep_view,
        size=5, alpha=0.75, color='slategrey')
dwell_ax.step(x='pooled_dwell_time', y='ecdf', source=pooled_source, view=pooled_view,
        legend='pooled data', line_width=2, color='dodgerblue')

# Plot the cutting idx
cut_ax.segment(x0=0, x1='value', y0='fate',  y1='fate', line_width=2, color='dodgerblue',
            source=fate_source, view=fate_view)
cut_ax.circle(x='value', y='fate', source=fate_source, view=fate_view, size=10, 
            fill_color='white', line_color='dodgerblue', line_width=2)

# Plot the cutting probabilities
pcut_ax.circle(x='p_cut', y='y', source=rep_pcut_source, color='slategrey',
                view=rep_pcut_view, size=4)

pcut_ax.segment(x0='low', x1='high', y0='y', y1='y', source=pcut_source, 
                color='dodgerblue', line_width=1, view=pcut_view)
pcut_ax.circle(x='p_cut', y='y', source=pcut_source, line_color='dodgerblue',
                view=pcut_view, fill_color='white', size=8, line_width=2)

# Plot the fractional looped time
floop_ax.circle(x='mean', y='y', source=floop_rep_source, view=floop_rep_view,
                color='slategrey', size=4)
floop_ax.circle(x='mean', y='y', source=floop_pooled_source, view=floop_view,
                fill_color='white', line_color='dodgerblue', line_width=2, size=10)

selector.js_on_change("value",dwell_cb)
dwell_ax.legend.location = 'bottom_right'
layout = bokeh.layouts.gridplot([[selector], [seq],
                                [dwell_ax],
                                [cut_ax, pcut_ax, floop_ax]])


# #############################################################################
#  THEME DETAILS
# #############################################################################

theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#f5e3b3',
                'outline_line_color': '#000000',
            },
            'Axis': {
            'axis_line_color': "black",
            'major_tick_out': 7,
            'major_tick_line_width': 0.75,
            'major_tick_line_color': "black",
            'minor_tick_line_color': "black",
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': None,
            },
            'Legend': {
                'background_fill_color': '#f5e3b3',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica',
                'offset': 2,
            }}}

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(layout)


