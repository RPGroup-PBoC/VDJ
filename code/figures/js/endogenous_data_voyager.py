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
from bokeh.models import ColumnDataSource, Div, LinearAxis, CustomJS, CDSView, Grid, GroupFilter, Band
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.embed import components
import vdj.io
import vdj.stats
import scipy.stats

# ##############################################################################
# DATA LOADING / TRANSFORMATION
# ##############################################################################
# Load the datasets
dwell_times = pd.read_csv('../../../data/compiled_dwell_times.csv')
dwell_times = dwell_times[dwell_times['mutant'] != 'analysis']
dwell_times = dwell_times[(~dwell_times['mutant'].str.startswith('12')) | (dwell_times['mutant']=='12SpacC1A')]
f_looped = pd.read_csv('../../../data/compiled_looping_fraction.csv')
f_looped = f_looped[f_looped['mutant'] != 'analysis']
f_looped = f_looped[(~f_looped['mutant'].str.startswith('12')) | (f_looped['mutant']=='12SpacC1A')]
fates = pd.read_csv('../../../data/compiled_cutting_events.csv')
fates = fates[fates['mutant'] != 'analysis']
fates = fates[(~fates['mutant'].str.startswith('12')) | (fates['mutant']=='12SpacC1A')]

# Load the sampling information
#samples = pd.read_csv('../../../data/pooled_model_samples.csv')
#samples = samples[(~samples['mutant'].str.startswith('12')) | (samples['mutant']=='12SpacC1A')]
#stats = pd.read_csv('../../../data/pooled_model_summary.csv')
#stats = stats[(~stats['mutant'].str.startswith('12')) | (stats['mutant']=='12SpacC1A')]

endog_names = {'WT12rss' : 'V4-57-1 (ref)',
                'DFL161' : 'DFL 16.1-5\'',
                'DFL1613' : 'DFL 16.1-3\'',
                '12SpacC1A' : 'V4-55'}

#dwell_times = dwell_times.replace({'mutant' : endog_names})
#f_looped = f_looped.replace({'mutant' : endog_names})
#fates = fates.replace({'mutant' : endog_names})
#samples = samples.replace({'mutant' : endog_names})
#stats = stats.replace({'mutant' : endog_names})

# Compute the statistics for the cutting probability and looping time
pcut_df = pd.DataFrame([])
floop_df = pd.DataFrame([])
for g, d in samples.groupby(['mutant']):
    ind = np.argmax(d['lp__'].values)
    # Compute p_loop
    _p_cut = d['r_cut'].values /\
                (d['r_cut'].values + d['k_unloop'].values)
    _floop = d['k_loop'].values /\
                (d['r_cut'].values + d['k_unloop'].values +\
                 d['k_loop'].values) 
    mean_pcut = np.mean(_p_cut)
    median_pcut = np.median(_p_cut) 
    mode_pcut = _p_cut[ind]
    hpd_min, hpd_max = vdj.stats.compute_hpd(_p_cut, 0.95)
    pcut_df = pcut_df.append({'mean':mean_pcut, 'median':median_pcut, 
                            'mode':mode_pcut, 'hpd_min':hpd_min, 
                            'hpd_max':hpd_max, 'mutant':g}, ignore_index=True)
    mean_floop = np.mean(_floop)
    median_floop = np.median(_floop) 
    mode_floop = _floop[ind]
    hpd_min, hpd_max = vdj.stats.compute_hpd(_floop, 0.95)
    floop_df = floop_df.append({'mean':mean_floop, 'median':median_floop, 
                            'mode':mode_floop, 'hpd_min':hpd_min, 
                            'hpd_max':hpd_max, 'mutant':g}, ignore_index=True)
pcut_df = pcut_df.replace({'mutant' : endog_names})
floop_df = floop_df.replace({'mutant' : endog_names})
# ##############################################################################
#  COMPUTING POOLED AND REPLICATE PROPERTIES 
# ##############################################################################
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
seq_df = seq_df.replace({'mutant' : endog_names})

# Compute the ECDFs for the dwell times
pooled_dwell_times= []
for g, d in dwell_times.groupby('mutant'):
    # Compute the pooled ecdf
    # x = np.sort(d['dwell_time_s'])
    # y = np.linspace(0, 1, len(d))
    hist, bins = np.histogram(d['dwell_time_min'], bins=100)
    hist = hist / hist.max()
    df = pd.DataFrame([])
    df['pooled_dwell_time'] = bins[:-1]
    df['ecdf'] = hist 
    df['mutant'] = g
    pooled_dwell_times.append(df)
pooled_dwell = pd.concat(pooled_dwell_times)
pooled_dwell = pooled_dwell.replace({'mutant' : endog_names})

# Compute the replicate ecdfs
rep_dwell_times = []
for g, d in dwell_times.groupby(['mutant', 'replicate']):
    # x  = np.sort(d['dwell_time_s'])
    # y = np.linspace(0, 1, len(d))
    hist, bins = np.histogram(d['dwell_time_min'], bins=100)
    hist = hist / hist.sum()
    df = pd.DataFrame([])
    df['dwell'] = bins[:-1]
    df['ecdf'] = hist 
    df['mutant'] = g[0] 
    df['replicate'] = g[1]
    rep_dwell_times.append(df)
rep_dwell = pd.concat(rep_dwell_times)
rep_dwell = rep_dwell.replace({'mutant' : endog_names})

# Compute the total bead counts for that particular mutant
pooled_fates = pd.DataFrame([])
for g, d in fates.groupby(['mutant']):
    pooled_fates = pooled_fates.append({'fate':'cut','value':d['n_cuts'].values.sum(),
                        'mutant':g}, ignore_index=True)
    pooled_fates = pooled_fates.append({'fate':'unlooped','value':d['n_beads'].values.sum() -\
                        d['n_cuts'].values.sum(), 'mutant':g},
                        ignore_index=True)
pooled_fates = pooled_fates.replace({'mutant' : endog_names})

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
cut_prob = cut_prob.replace({'mutant' : endog_names})


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
rep_cut_prob = rep_cut_prob.replace({'mutant' : endog_names})

# Compute the fractional time in looped state for all data.
floop_pooled = f_looped.groupby(['mutant'])['fraction_looped'].agg(('mean', 'sem')).reset_index()
floop_pooled['y'] = 0
floop_pooled = floop_pooled.replace({'mutant' : endog_names})
floop_rep = f_looped.groupby(['mutant', 'replicate'])['fraction_looped'].agg(('mean', 'sem')).reset_index()
floop_rep['y'] = np.random.normal(0, 0.01, len(floop_rep))
floop_rep = floop_rep.replace({'mutant' : endog_names})

# %%
# ##############################################################################
# FIGURE SOURCE AND VIEW DEFINITIONS
# ##############################################################################
bokeh.plotting.output_file('./data_voyager.html')

# Instantiate the dropdown menu
selector = Select(title='Mutant', value='V4-57-1 (ref)', 
                  options=list(floop_rep.mutant.unique()))

# Define the filter for mutant selection
filter = GroupFilter(column_name="mutant", group="V4-57-1 (ref)")


##################### Generate the sources #####################################

# Dwell times
rep_source = ColumnDataSource(rep_dwell)
pooled_source = ColumnDataSource(pooled_dwell)

# Modeled dwell times
#dwell_cdf_source = ColumnDataSource(dwell_cdfs)

# Bead fate
fate_source = ColumnDataSource(pooled_fates)

# Cutting probabilities
pcut_source = ColumnDataSource(cut_prob)
rep_pcut_source = ColumnDataSource(rep_cut_prob)

# Modeled cutting probabily
#pcut_model_source = ColumnDataSource(pcut_df)

# Looping fraction
floop_pooled_source = ColumnDataSource(floop_pooled)
floop_rep_source = ColumnDataSource(floop_rep)

# Modeled looping fraction
#floop_model_source = ColumnDataSource(floop_df)

# Sequences
seq_source = ColumnDataSource(seq_df)


######################## Assign Data Views #####################################

# Dwell times
pooled_view = CDSView(source=pooled_source, filters=[filter])
rep_view = CDSView(source=rep_source, filters=[filter])

# Modelled dwell times
#dwell_cdf_view = CDSView(source=dwell_cdf_source, filters=[filter])

# Bead Fates
fate_view = CDSView(source=fate_source, filters=[filter])

# Cutting probabilities
pcut_view = CDSView(source=pcut_source, filters=[filter])
rep_pcut_view = CDSView(source=rep_pcut_source, filters=[filter])

# Modeled p_cut view
pcut_model_view = CDSView(source=pcut_model_source, filters=[filter])

# Looping Fraction
floop_view = CDSView(source=floop_pooled_source, filters=[filter])
floop_rep_view = CDSView(source=floop_rep_source, filters=[filter])

# Modeled floop view
floop_model_view = CDSView(source=floop_model_source, filters=[filter])

# Sequence display
seq_view = CDSView(source=seq_source, filters=[filter])

# ##############################################################################
#  Js CALLBACK DEFINITION
# ##############################################################################

# Define the callback args and callback code
cb_args = {'pv':pooled_view, 'rv':rep_view, 'fv':fate_view, 'pcv':pcut_view, 
           'rpcv':rep_pcut_view, 'flv':floop_view, 'flrv':floop_rep_view,
           'seqv':seq_view, 'cdfv':dwell_cdf_view, 'cflv':floop_model_view,
           'cpcv': pcut_model_view,
          'sel':selector, 'filter':filter,
          'ps':pooled_source, 'rs':rep_source, 'fs':fate_source, 
          'pcs':pcut_source, 'rpcs':rep_pcut_source, 'fls':floop_pooled_source,
          'flrs':floop_rep_source, 'seqs':seq_source, 'cdfs':dwell_cdf_source,
          'cfls': floop_model_source, 'cpcs':pcut_model_source}
dwell_cb = CustomJS(args=cb_args, 
    code="""
    var mut = sel.value
    filter.group = mut;

    // Define the filters
    pv.filters[0] = filter; // Pooled dwell times
    rv.filters[0] = filter; //  Replicate dwell times
    seqv.filters[0] = filter; // Sequences
    fv.filters[0] = filter;  // Bead fates
    pcv.filters[0] = filter; // Pooled cutting probability 
    rpcv.filters[0] = filter; // Replicate cutting probabily
    flv.filters[0] = filter; // Pooled looping fraction
    flrv.filters[0] = filter; // replicate looping fraction
    cdfv.filters[0] = filter; // Modeled dwell time
    cflv.filters[0] = filter; // Modeled looping fraction
    cpcv.filters[0] = filter; // Modeled cutting probability
    
    // Update the view on the data
    seqs.data.view = seqs; // Sequence source
    fls.data.view = fv; // Pooled looping fraction source 
    flrs.data.view = fv; // Replicate looping fraction source
    ps.data.view = pv; // Pooled dwell time source
    rs.data.view = rv; // Replicate dwell time source
    pcs.data.view = pv; // Pooled cutting probability source
    rpcs.data.view = pv; // Replicate cutting probability source
    fs.data.view = fv; // Bead fate source
    cdfs.data.view = cdfv; // Modeled dwell times
    cfls.data.view = cflv; // Modeled looping fraction
    cpcs.data.view = cpcv; // modeled cutting probability

    // Push changes to the plot 
    rs.change.emit(); // Replicate dwell times
    ps.change.emit(); // Pooled dwell times
    pcs.change.emit(); // Pooled cutting probabilities
    rpcs.change.emit(); // Replicate cutting probabilities
    fs.change.emit(); // Bead fates
    fls.change.emit(); // Pooled looping fraction
    flrs.change.emit(); // Replicate looping fraction
    seqs.change.emit(); // Sequence display
    cdfs.change.emit(); // Modeled dwell times
    cfls.change.emit(); // Modeled looping fraction
    cpcs.change.emit(); // Modeled cutting probability
""")

# ##############################################################################
# DEFINE THE FIGURE CANVAS
# ##############################################################################
# Dwell time figure canvas
dwell_ax = bokeh.plotting.figure(width=800, height=400, 
                           x_axis_label='dwell time [min]',
                           y_axis_label='cumulative distribution',
                           title=f'paired complex dwell time',
                           x_range=[np.min(dwell_times['dwell_time_min']), 
                                    np.max(dwell_times['dwell_time_min'])],
                           x_axis_type='log')

_fates = ['unlooped', 'cut']

# Cutting fraction canvas
cut_ax = bokeh.plotting.figure(width=200, height=150, 
                               x_axis_label='number of events',
                               y_range=_fates, title='bead fate')
                               
# Cutting Probability Canvas
pcut_ax = bokeh.plotting.figure(width=300, height=150, x_range=[-0.1, 1.1],
                               x_axis_label='probability', y_range=[-0.1, 0.1],
                title=f'cutting probability')
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

# ##############################################################################
# PLOTTING
# ##############################################################################
# Plot the dwell time distributions
# dwell_ax.circle(x=[], y=[], line_width=2, color='slategrey', legend='replicate')
# dwell_ax.circle(x='dwell', y='ecdf', source=rep_source, view=rep_view,
#         size=5, alpha=0.75, color='slategrey')
dwell_ax.step(x='pooled_dwell_time', y='ecdf', source=pooled_source, view=pooled_view,
        legend='pooled data', line_width=2, color='dodgerblue')
    
dwell_ax.quad(top='hpd_max', bottom='hpd_min', left='left', 
            right='right', color='tomato', alpha=0.5,
            source=dwell_cdf_source,
            view=dwell_cdf_view, legend='fit')

# Plot the cutting idx
cut_ax.segment(x0=0, x1='value', y0='fate',  y1='fate', line_width=2, color='dodgerblue',
            source=fate_source, view=fate_view)
cut_ax.circle(x='value', y='fate', source=fate_source, view=fate_view, size=10, 
            fill_color='white', line_color='dodgerblue', line_width=2)

# Plot the cutting probabilities
pcut_ax.circle(x='p_cut', y='y', source=rep_pcut_source, color='slategrey',
                view=rep_pcut_view, size=4)
pcut_ax.circle(x='p_cut', y=0.05, source=pcut_source, line_color='dodgerblue',
                view=pcut_view, fill_color='white', size=8, line_width=2)

# Plot the modeled cutting probabilities
pcut_ax.square(x='median', y=-0.05, source=pcut_model_source, view=pcut_model_view,
               fill_color='white', size=8, line_color='tomato')
pcut_ax.segment(x0='hpd_min', x1='hpd_max', y0=-0.05, y1=-0.05,
               source=pcut_model_source, line_width=2, color='tomato',
               view=pcut_model_view)

# Plot the fractional looped time
floop_ax.circle(x='mean', y='y', source=floop_rep_source, view=floop_rep_view,
                color='slategrey', size=4)
floop_ax.circle(x='mean', y=0.05, source=floop_pooled_source, view=floop_view,
                fill_color='white', line_color='dodgerblue', line_width=2,
                size=10)
                
# Modeled looping fraction
floop_ax.square(x='median', y=-0.05, source=floop_model_source, view=floop_model_view,
               fill_color='white', size=8, line_color='tomato')
floop_ax.segment(x0='hpd_min', x1='hpd_max', y0=-0.05, y1=-0.05,
               source=floop_model_source, line_width=2, color='tomato',
               view=floop_model_view)


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



