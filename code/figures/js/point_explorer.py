#%%
# -*- coding: utf-8 -*- #%%
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, 
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool,
                          LinearColorMapper, TapTool, RadioButtonGroup)
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select 
from bokeh.embed import components
import vdj.io
import vdj.stats
from bokeh.transform import transform
import bokeh.palettes
# bokeh.io.output_notebook()

# Load the necessary data sets
dwell_times = pd.read_csv('../../../data/compiled_dwell_times.csv')
posteriors = pd.read_csv('../../../data/pooled_cutting_probability_posteriors.csv')
loops = pd.read_csv('../../../data/compiled_looping_events.csv')
pcuts = pd.read_csv('../../../data/pooled_cutting_probability.csv')

# Keep the data from the reference sequence
dwell_ref = dwell_times[(dwell_times['mutant']=='WT12rss') & 
                       (dwell_times['hmgb1']==80) & (dwell_times['salt']=='Mg')]
cut_ref = dwell_ref[dwell_ref['cut']==1]
post_ref = posteriors[(posteriors['mutant']=='WT12rss') & 
                             (posteriors['hmgb1']==80) & (posteriors['salt']=='Mg')]
loops_ref = loops[(loops['mutant']=='WT12rss') & 
                    (loops['hmgb1']==80) & (loops['salt']=='Mg')]
pcut_ref = pcuts[(pcuts['mutant']=='WT12rss') &
            (pcuts['hmgb1']==80) &(pcuts['salt']=='Mg')]

# Identify the endogenous muts
dfs  = []
for i, df in enumerate([dwell_times, posteriors, loops, pcuts]):
    for g, d in df.groupby(['mutant']):
        n_muts = vdj.io.mutation_parser(g)['n_muts']
        if n_muts > 1: 
            mut_class = 'endogenous'
        else:
            mut_class = 'point' 
       
        df.loc[df['mutant']==g, 'class'] = mut_class
    df = df[(df['class']=='point') & (df['salt']=='Mg') & (df['hmgb1']==80)]
    dfs.append(df)
dwell_times, posteriors, loops, pcuts = dfs
cut_dwells = dwell_times[dwell_times['cut']==1]
post_dist_source = ColumnDataSource(posteriors)
#%%
# ##############################################################################
# COMPUTE POOLED v REPLICATE LOOP FREQUENCIES
# ##############################################################################
pooled_df = pd.DataFrame()
rep_df = pd.DataFrame()
for g, d in loops.groupby('mutant'):
    pooled_df = pooled_df.append({'mutant': g, 
        'loops_per_bead':d['n_loops'].sum() / len(d['bead_idx']),
        'n_beads':len(d['bead_idx']), 'n_loops':d['n_loops'].sum()},
        ignore_index=True)
    for _g, _d in d.groupby(['replicate']):
        rep_df = rep_df.append({'mutant':g,
            'loops_per_bead':_d['n_loops'].sum() / len(_d['bead_idx']),
            'n_beads':len(_d['bead_idx']), 'n_loops':_d['n_loops'].sum()},
            ignore_index=True)
pooled_df['y'] =0.5 
rep_df['y'] = np.random.normal(0.5,0.05, len(rep_df))

# Assemble into sources
pooled_dist_source = ColumnDataSource(pooled_df)
rep_dist_source = ColumnDataSource(rep_df)

# Do the same for the reference
pooled_ref = pd.DataFrame()#
rep_ref = pd.DataFrame()
for g, d in loops_ref.groupby('mutant'):
    pooled_ref = pooled_ref.append({'mutant': g, 
        'loops_per_bead':d['n_loops'].sum() / len(d['bead_idx']),
        'n_beads':len(d['bead_idx']), 'n_loops':d['n_loops'].sum()},
        ignore_index=True)
    for _g, _d in d.groupby(['replicate']):
        rep_ref = rep_ref.append({'mutant':g,
            'loops_per_bead':_d['n_loops'].sum() / len(_d['bead_idx']),
            'n_beads':len(_d['bead_idx']), 'n_loops':_d['n_loops'].sum()},
            ignore_index=True)
pooled_ref['y'] = 0
rep_ref['y'] = np.random.normal(0,0.05, len(rep_ref))
#%%
# ##############################################################################
# GENERATE HISTOGRAMS OF DWELL TIMES 
# ##############################################################################
# Generate the histogrammed dwell times.
bins = np.linspace(0, dwell_times['dwell_time_min'].max(), 25)
dfs = []
for source in [dwell_times, cut_dwells]: 
    bin_dfs = []
    for g, d in source.groupby('mutant'):
        hist, bins = np.histogram(d['dwell_time_min'], bins=bins)
        _df = pd.DataFrame()
        _df['top'] = hist
        _df['bottom'] = 0
        _df['left'] = bins[1:]
        _df['right'] = bins[:-1]
        _df['mutant'] = g
        bin_dfs.append(_df)
    dwell_hist = pd.concat(bin_dfs)
    dfs.append(dwell_hist)
dwell_hist, cut_hist = dfs

# Assemble into sources
dwell_hist_source = ColumnDataSource(dwell_hist)
cut_hist_source = ColumnDataSource(cut_hist)

# Do the same for the reference sequence
dfs = []
for source in [dwell_ref, cut_ref]: 
    bin_dfs = []
    for g, d in source.groupby('mutant'):
        hist, bins = np.histogram(d['dwell_time_min'], bins=bins)
        _df = pd.DataFrame()
        _df['top'] = hist
        _df['bottom'] = 0
        _df['left'] = bins[1:]
        _df['right'] = bins[:-1]
        _df['mutant'] = g
        bin_dfs.append(_df)
    dwell_hist = pd.concat(bin_dfs)
    dfs.append(dwell_hist)
dwell_hist_ref, cut_hist_ref = dfs

#%%
# ##############################################################################
# COMPUTE MEDIAN DWELL TIME DATA
# ##############################################################################
median_dwell_ref = dwell_ref['dwell_time_min'].median()
# Set up a test sequence map
dwell_mat = pd.DataFrame()
ref_seq = vdj.io.endogenous_seqs()['WT12rss'][1]
nt_idx = vdj.io.nucleotide_idx()
# def parse_mutation(0)
for g, d in dwell_times.groupby(['mutant']):
    mut_seq = vdj.io.mutation_parser(g)
    if g != 'WT12rss':
        # Parse the mutation
        loc = np.where(mut_seq['seq_idx']!=ref_seq)[0][0]
        base = mut_seq['seq'][loc]

        # Compute the median dwell time
        med_dwell = d['dwell_time_min'].median()
        diff = med_dwell - median_dwell_ref
        dwell_mat = dwell_mat.append(
                {'mutant': g,
                'pos': loc,
                'base_idx': nt_idx[base],
                'base': base,
                'diff': diff,
                'med_dwell':med_dwell}, ignore_index=True)
dwell_source = ColumnDataSource(dwell_mat)

#%%
# ##############################################################################
# COMPUTE THE DIFFERENCE IN LOOP FREQUENCY
# ##############################################################################
pooled_loop_mat = pd.DataFrame()
rep_loop_mat = {}
for g, d in pooled_df.groupby(['mutant']):
    mut_seq = vdj.io.mutation_parser(g)
    if g != 'WT12rss':
        # Parse the mutation
        loc = np.where(mut_seq['seq_idx']!=ref_seq)[0][0]
        base = mut_seq['seq'][loc]

        # Populate the data frame
        diff = d['loops_per_bead'].values[0] - pooled_ref ['loops_per_bead'].values[0]
        pooled_loop_mat= pooled_loop_mat.append(
                 {'mutant': g,
                 'pos': loc,
                 'base_idx': nt_idx[base],
                 'loops_per_bead': d['loops_per_bead'].values[0],
                 'base': base,
                 'diff': diff}, ignore_index=True)
loop_source = ColumnDataSource(pooled_loop_mat)
#%%
# ##############################################################################
# COMPUTE THE DIFFERENCE IN CUTTING PROBABILITY
# ##############################################################################
mean_cut_ref = pcut_ref['mode'].values[0]
pcut_mat = pd.DataFrame()
for g, d in pcuts.groupby(['mutant']):
    mut_seq = vdj.io.mutation_parser(g)
    if g != 'WT12rss':
        # Parse the mutation
        loc = np.where(mut_seq['seq_idx']!=ref_seq)[0][0]
        base = mut_seq['seq'][loc]

        # Populate the data frame
        val = d['mode'].values[0]
        diff = val - mean_cut_ref
        pcut_mat = pcut_mat.append(
                 {'mutant': g,
                 'pos': loc,
                 'base_idx': nt_idx[base],
                 'pcut': val,
                 'base': base,
                 'diff': diff}, ignore_index=True)
pcut_source = ColumnDataSource(pcut_mat)
# %%
# Set up the matrices
ax_loop_mat = bokeh.plotting.figure(height=120, x_range=[-1, 28], tools=[''])
ax_loop = bokeh.plotting.figure(height=120, x_axis_label='paired complexes per bead',
            y_range=[-0.3, 0.8])
                               
ax_dwell_mat = bokeh.plotting.figure(height=120, x_range=[-1, 28], tools=[''])

ax_dwell = bokeh.plotting.figure(height=200, 
        x_axis_label='paired-complex dwell time [min]', y_axis_label='number of observations',
        tools=[''])
ax_cut_mat = bokeh.plotting.figure(height=120, x_range=[-1, 28], tools=[''])
                           

ax_cut = bokeh.plotting.figure(height=200, 
    x_axis_label='cleavage probability', y_axis_label='posterior probability',
    tools=[''])


# Insert interactivity
mut_filter = GroupFilter(column_name="mutant", group='')
rep_view = CDSView(source=rep_dist_source, filters=[mut_filter])
pooled_view = CDSView(source=pooled_dist_source, filters=[mut_filter])
dwell_view = CDSView(source=dwell_hist_source, filters=[mut_filter])
cut_view = CDSView(source=cut_hist_source, filters=[mut_filter])
post_view = CDSView(source=post_dist_source, filters=[mut_filter])



# Set the ticks to the reference sequence
tooltips = [('mutation', '@mutant'), ('difference from ref.', '@diff')]
hover_cb = CustomJS(args={'mut_filter':mut_filter, 'mut_source':loop_source,
                          'rep_view':rep_view, 'pooled_view':pooled_view,
                          'dwell_view':dwell_view, 'cut_view':cut_view,
                          'post_view':post_view, 'rep_source':rep_dist_source,
                          'pooled_source':pooled_dist_source, 
                          'dwell_source': dwell_hist_source,
                          'cut_source': cut_hist_source,
                          'post_source':post_dist_source}, 
code="""
var mut_ind = cb_data.index['1d'].indices[0];
var mut = mut_source.data['mutant'][mut_ind];
mut_filter.group = mut;
console.log(mut_filter.group)
var views = [rep_view, pooled_view, dwell_view, cut_view, post_view];
var sources = [rep_source, pooled_source, dwell_source, cut_source, post_source];
for (var i = 0; i < views.length; i++) {
    views[i].filters[0] = mut_filter;
    sources[i].data.view = views[i];
    sources[i].change.emit();
}
""")

for a, s  in zip([ax_loop_mat, ax_dwell_mat, ax_cut_mat], [loop_source, dwell_source, pcut_source]):
    hover = HoverTool(renderers=a, tooltips=tooltips, callback=hover_cb)
    a.add_tools(hover)
    a.yaxis.ticker = [0, 1, 2, 3]
    a.xaxis.ticker = np.arange(0, 29, 1)
    ylab = {int(i):nt_idx[i] for i in range(4)}
    xlab = {int(i):b for i, b in zip(np.arange(0, 27, 1), str(vdj.io.endogenous_seqs()['WT12rss'][0]))}
    a.yaxis.major_label_overrides = ylab
    a.xaxis.major_label_overrides = xlab

ax_loop.yaxis.visible = False


# Define titles
ax_loop_mat.title.text = "paired complex formation frequency"
ax_dwell_mat.title.text = "paired complex dwell time"
ax_cut_mat.title.text = "paired complex cleavage probability"

# Define the layout
# Define whether hover or selection should be used. 
button = RadioButtonGroup(labels=["Display on Hover", "Display on Click"], 
            active=0)
lay = bokeh.layouts.column(ax_loop_mat, ax_loop, ax_dwell_mat, 
                            ax_dwell, ax_cut_mat, ax_cut)

# Define the color palettes
loop_palette = bokeh.palettes.RdBu11
dwell_palette = bokeh.palettes.RdGy11
cut_palette = bokeh.palettes.Spectral11
loop_color = LinearColorMapper(palette=loop_palette, low=-0.2, high=0.2)
dwell_color = LinearColorMapper(palette=dwell_palette, low=-2, high=2)
cut_color = LinearColorMapper(palette=cut_palette, low=-0.5, high=0.6)


# Populate the matrices.
loop_fig = ax_loop_mat.rect('pos', 'base_idx', width=1, height=1, source=loop_source, 
        fill_color=transform('diff', loop_color))
dwell_fig = ax_dwell_mat.rect('pos', 'base_idx', width=1, height=1, source=dwell_source, 
        fill_color=transform('diff', dwell_color))
cut_fig = ax_cut_mat.rect('pos', 'base_idx', width=1, height=1, source=pcut_source, 
        fill_color=transform('diff', cut_color))

# Add the reference features
ax_loop.triangle('loops_per_bead', 'y', color='grey', alpha=0.5, source=rep_ref, 
                   size=8)
ax_loop.circle('loops_per_bead', 'y', line_color='grey', fill_color='white', 
                alpha=0.5, source=pooled_ref, size=10, line_width=2)

ax_dwell.quad(bottom=0, left='left', right='right', top='top', color='grey', 
            alpha=0.3, source=dwell_hist_ref)
ax_dwell.quad(bottom=0, left='left', right='right', top='top', color=None, 
            hatch_alpha=0.5, source=cut_hist_ref, hatch_pattern='/', hatch_color='grey',
            line_color='grey', line_alpha=0.3)

ax_cut.line('probability', 'posterior', source=post_ref, color='grey',
            alpha=0.3)
ax_cut.varea('probability', 0, 'posterior', source=post_ref, fill_color='grey',
            alpha=0.3)


# Add the point mutant features.
ax_loop.triangle('loops_per_bead', 'y', color='dodgerblue', alpha=0.5, 
                 source=rep_dist_source, view=rep_view)
ax_loop.circle('loops_per_bead', 'y', line_color='dodgerblue', alpha=0.75,
                fill_color='white', source=pooled_dist_source, view=pooled_view)
ax_dwell.quad(bottom=0, left='left', right='right', top='top',
                 color='dodgerblue', alpha=0.5, source=dwell_hist_source,
                 view=dwell_view)
ax_dwell.quad(bottom=0, left='left', right='right', top='top',
                 line_color='navy', hatch_color='navy', hatch_pattern='\\',
                 alpha=0.5, source=cut_hist_source,
                 view=cut_view)
ax_cut.line('probability', 'posterior', source=post_dist_source, color='dodgerblue',
            view=post_view)
ax_cut.varea('probability', y1=0, y2='posterior', source=post_dist_source, 
            fill_color='dodgerblue', view=post_view, alpha=0.5)


for a, s  in zip([ax_loop_mat, ax_dwell_mat, ax_cut_mat], [loop_fig, dwell_fig, cut_fig]):
    hover = HoverTool(renderers=[s], tooltips=tooltips, callback=hover_cb)
    a.add_tools(hover)

# Plot x's for virgin mutation
x, y = [], []
for p in range(28):
   for b in range(4):
      if (len(pcut_mat[(pcut_mat['pos'] == p) &
         (pcut_mat['base'] == nt_idx[b])]) == 0) & (ref_seq[p] != b):
          x.append(p)
          y.append(b)

for a in [ax_loop_mat, ax_dwell_mat, ax_cut_mat]:
    a.x(x, y, color='grey')
    # Fill in the wild-type positions
    a.rect(np.arange(0, 29, 1), ref_seq, width=1, height=1, fill_color='white',
            line_color='slategrey')
    a.circle(np.arange(0, 29, 1), ref_seq,  fill_color='#f5e3b3', alpha=0.75,
              line_color='slategrey', line_width=1, size=10)



# ##############################################################################
# INTERACTIVITY DEFINITION
# ##############################################################################
# Define the callback for tap 
# Apply the theme. 
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
                'text_font_style': 'bold',
                'align': 'center',
                'text_font': 'Helvetica',

                'offset': 2,
            }}}
bokeh.plotting.output_file('./point_voyager.html')

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(lay)


#%%
