#%%
# -*- coding: utf-8 -*- #%%
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.events import Tap
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, 
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool,
                          LinearColorMapper, TapTool, RadioButtonGroup,
                          ColorBar, FixedTicker, Button, Segment)
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
loops = pd.read_csv('../../../data/compiled_loop_freq_bs.csv')
pcuts = pd.read_csv('../../../data/pooled_cutting_probability.csv')
pcuts.rename(columns={'n_beads': 'n_loops'}, inplace=True)
# Keep the data from the reference sequence
dwell_ref = dwell_times[(dwell_times['mutant']=='WT12rss') & 
                       (dwell_times['hmgb1']==80) & (dwell_times['salt']=='Mg')]
cut_ref = dwell_ref[dwell_ref['cut']==1]
unlooped_ref = dwell_ref[dwell_ref['cut']==0]
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
        if 'cod' in g.lower():
            mut_class = 'flank'
        elif (n_muts > 1) | (g == 'V4-55'): 
            mut_class = 'endogenous'
        else:
            mut_class = 'point' 
        df.loc[df['mutant']==g, 'class'] = mut_class
    df = df[(df['class']=='point') & (df['salt']=='Mg') & (df['hmgb1']==80)]
    dfs.append(df)
dwell_times, posteriors, loops, pcuts = dfs
cut_dwells = dwell_times[dwell_times['cut']==1]
unlooped_dwells = dwell_times[dwell_times['cut']==0]
post_dist_source = ColumnDataSource(posteriors)
#%%
loop_source = ColumnDataSource(loops)
loop_ref_source = ColumnDataSource(loops_ref)

#%%
# ##############################################################################
# GENERATE ECDFS OF DWELL TIMES 
# ##############################################################################

# Generate the histogrammed dwell times.
# bins = np.linspace(0, dwell_times['dwell_time_min'].max(), 25)
dfs = []
for source in [dwell_times, cut_dwells, unlooped_dwells]: 
    bin_dfs = []
    for g, d in source.groupby('mutant'):
        x, y = np.sort(d['dwell_time_min'].values), np.arange(0, len(d), 1) / len(d)
        y[-1] = 1
        _df = pd.DataFrame()
        _df['dwell'] = x
        _df['ecdf'] = y
        _df['mutant'] = g
        bin_dfs.append(_df)
    dwell_dist = pd.concat(bin_dfs)
    dfs.append(dwell_dist)
dwell_dist, cut_dist, unlooped_dist = dfs

# Assemble into sources
dwell_dist_source = ColumnDataSource(dwell_dist)
cut_dist_source = ColumnDataSource(cut_dist)
unlooped_dist_source = ColumnDataSource(unlooped_dist)

# Do the same for the reference sequence
dfs = []
for source in [dwell_ref, cut_ref, unlooped_ref]: 
    bin_dfs = []
    for g, d in source.groupby('mutant'):
        x, y = np.sort(d['dwell_time_min'].values), np.arange(0, len(d), 1) / len(d) 
        y[-1] = 1
        _df = pd.DataFrame()
        _df['dwell'] = x
        _df['ecdf'] = y
        _df['mutant'] = g
        bin_dfs.append(_df)
    dwell_dist = pd.concat(bin_dfs)
    dfs.append(dwell_dist)
dwell_dist_ref, cut_dist_ref, unlooped_dist_ref = dfs

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
    if (g != 'WT12rss') & ('cod' not in g.lower()):
        # Parse the mutation
        loc = np.where(mut_seq['seq_idx']!=ref_seq)[0][0]
        try:
            base = mut_seq['seq'][loc]
        except TypeError:
            base = mut_seq['seq']

        # Compute the median dwell time
        med_dwell = d['dwell_time_min'].median()
        diff = med_dwell - median_dwell_ref
        dwell_mat = dwell_mat.append(
                {'mutant': g,
                'pos': loc,
                'base_idx': nt_idx[base],
                'base': base,
                'diff': diff,
                'med_dwell':med_dwell,
                'n_loops': int(len(d))}, ignore_index=True)
dwell_source = ColumnDataSource(dwell_mat)

#%%
# ##############################################################################
# COMPUTE THE DIFFERENCE IN LOOP FREQUENCY
# ##############################################################################
pooled_loop_mat = pd.DataFrame()
for g, d in loops.groupby(['mutant']):
    mut_seq = vdj.io.mutation_parser(g)
    if (g != 'WT12rss') & ('cod' not in g.lower()):
        # Parse the mutation
        loc = np.where(mut_seq['seq_idx']!=ref_seq)[0][0]
        base = mut_seq['seq'][loc]

        # Populate the data frame
        diff = d['loops_per_bead'].values[0] - loops_ref ['loops_per_bead'].values[0]
        pooled_loop_mat= pooled_loop_mat.append(
                 {'mutant': g,
                 'pos': loc,
                 'base_idx': nt_idx[base],
                 'loops_per_bead': d['loops_per_bead'].values[0],
                 'base': base,
                 'diff': diff,
                 'n_beads': int(d['n_beads'].unique()),
                 'n_loops': int(d['n_loops'].unique())}, ignore_index=True)
loop_source = ColumnDataSource(pooled_loop_mat)
#%%
# ##############################################################################
# COMPUTE THE DIFFERENCE IN CUTTING PROBABILITY
# ##############################################################################
mean_cut_ref = pcut_ref['mode'].values[0]
pcut_mat = pd.DataFrame()
for g, d in pcuts.groupby(['mutant']):
    mut_seq = vdj.io.mutation_parser(g)
    if (g != 'WT12rss') & ('cod' not in g.lower()):
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
                 'diff': diff,
                 'n_loops': int(d['n_loops'].unique()),
                 'n_cuts': int(d['n_cuts'].unique())}, ignore_index=True)
pcut_source = ColumnDataSource(pcut_mat)
# %%
# Set up the matrices
ax_loop_mat = bokeh.plotting.figure(height=120, width=600, x_range=[-1, 28], tools=['tap'],
            toolbar_location=None)
ax_loop = bokeh.plotting.figure(height=120, width=600, x_axis_label='paired complexes per bead',
            y_range=[-0.8, 0.8], x_range=[-0.1, 0.95], toolbar_location=None)

ax_dwell_mat = bokeh.plotting.figure(height=120, width=600, x_range=[-1, 28], tools=['tap'],
                toolbar_location=None)

ax_dwell_unlooped = bokeh.plotting.figure(height=200, width=185,
        x_axis_label='dwell time [min]', y_axis_label='ECDF',
        tools=[''], toolbar_location=None, x_axis_type='log', title='unlooped PCs',
        x_range=[0.50, 80])
ax_dwell_cut = bokeh.plotting.figure(height=200, width=185, 
        x_axis_label='dwell time [min]', y_axis_label='ECDF',
        tools=[''], toolbar_location=None, x_axis_type='log', title='cleaved PCs',
        x_range=[0.50, 80])

ax_dwell_all = bokeh.plotting.figure(height=200, width=185,
        x_axis_label='dwell time [min]', y_axis_label='ECDF',
        tools=[''], toolbar_location=None, x_axis_type='log', title='all PCs',
        x_range=[0.50, 80])
ax_cut_mat = bokeh.plotting.figure(height=120, width=600, x_range=[-1, 28], tools=['tap'],
            toolbar_location=None)
ax_cut = bokeh.plotting.figure(height=200, width=600, 
    x_axis_label='cleavage probability', y_axis_label='posterior probability',
    tools=[''], toolbar_location=None)



# Add a blank legend plot. 
ax_leg = bokeh.plotting.figure(height=50, width=600, tools=[''], toolbar_location=None)
ax_leg2 = bokeh.plotting.figure(height=60, width=600, tools=[''], toolbar_location=None)
ax_leg.rect([], [], width=1, height=1, fill_color='white', line_color='black', legend='reference nucleotide')
ax_leg.circle([], [], fill_color='#f5e3b3', line_color='black', size=15, legend='reference nucleotide')
ax_leg.x([], [], color='black', size=10, legend='not measured')
ax_leg.line([], [], color='slategrey', legend='reference data')
ax_leg.circle([], [], color='slategrey', legend='reference data')
ax_leg.triangle([], [], color='slategrey', legend='reference data')
ax_leg.title.text_font_style = "normal"
ax_leg.legend.spacing = 50 
ax_leg.legend.location='center'
ax_leg.legend.orientation = 'horizontal'
ax_leg.legend.background_fill_color='white'

for a in [ax_leg, ax_leg2]:
    a.outline_line_color = None
    a.background_fill_color = 'white'
    a.xaxis.visible = False
    a.yaxis.visible = False
    a.xgrid.visible = False
    a.ygrid.visible = False

# Insert interactivity
mut_filter = GroupFilter(column_name="mutant", group='')
loop_view = CDSView(source=loop_source, filters=[mut_filter])
dwell_view = CDSView(source=dwell_dist_source, filters=[mut_filter])
unlooped_view = CDSView(source=unlooped_dist_source, filters=[mut_filter])
cut_view = CDSView(source=cut_dist_source, filters=[mut_filter])
post_view = CDSView(source=post_dist_source, filters=[mut_filter])

colors1 = bokeh.palettes.Blues9[1:-2]
colors2 = bokeh.palettes.Greys9[1:-2]


perc_source = []
perc_view = []
percs = list(np.sort(loops['percentile'].unique()))
percs.reverse()
ax_loop.triangle(x='loops_per_bead', y=-0.5, source=loop_ref_source,
                fill_color='white', line_color='slategrey', 
                size=10, level='overlay', legend='observed frequency')
ax_loop.triangle(x='loops_per_bead', y=0.5, source=loop_source,
                view=loop_view, fill_color='white', line_color='dodgerblue', 
                size=10, level='overlay', legend='observed frequency')

for i, p in enumerate(percs):
    d = loops[loops['percentile']==p]
    d_ref = ColumnDataSource(loops_ref[loops_ref['percentile']==p])
    _source = ColumnDataSource(d)
    _view = CDSView(source=_source, filters=[mut_filter])
    perc_source.append(_source)
    perc_view.append(_view)

    band = Segment(x0='low', x1='high', y0=0.5, y1=0.5, 
                    line_color=colors1[-1 * (i+1)], line_width=25) 
    ref_band = Segment(x0='low', x1='high', y0=-0.5, y1=-0.5, 
                    line_color=colors2[-1 * (i+1)], line_width=25) 
    
    ax_loop.add_glyph(_source, band, view=_view)
    ax_loop.add_glyph(d_ref, ref_band)




# Set the ticks to the reference sequence

loop_sel_code = """
var mut_ind = loop_source.selected['1d'].indices[0];
var mut = loop_source.data['mutant'][mut_ind];
"""

dwell_sel_code = """
var mut_ind = dwell_source.selected['1d'].indices[0];
var mut = dwell_source.data['mutant'][mut_ind];
"""

cut_sel_code = """
var mut_ind = cut_source.selected['1d'].indices[0];
var mut = cut_source.data['mutant'][mut_ind];
"""

sel_code = """
var sources = [loop_source, cut_source, dwell_source];
for (var i = 0; i < sources.length; i++) {
    sources[i].selected['1d'].indices[0] = mut_ind;
    sources[i].change.emit();
} 
"""

reset_code = """
    var mut = '';
    var plots = [loop_mat, dwell_mat, cut_mat];
    for (var i = 0; i < plots; i++) {
        plots.reset.emit();
       }
     """

draw = """
mut_filter.group = mut;

// Update the percentiles
for (var i = 0; i < percentile_source.length; i++ ) {
    var perc = percentile_source[i];
    percentile_view[i].filters[0] = mut_filter
    perc.data.view = percentile_view;
    perc.change.emit();
}

views = [loop_view, dwell_view, cut_view, unlooped_view, pcut_view];
data = [loop_source, dwell_data, cut_data, unlooped_data, pcut_data];

for (var i = 0; i < views.length; i++) {
    views[i].filters[0] = mut_filter;
    data[i].data.view = views[i];
    data[i].change.emit();
}
"""

args = {'mut_filter':mut_filter, 
        'loop_source':loop_source, 
        'dwell_source':dwell_source,
        'cut_source':pcut_source, 
        'loop_view':loop_view, 
        'dwell_view':dwell_view,
        'cut_view':cut_view,
        'pcut_view':post_view,
        'unlooped_view':unlooped_view,
        'dwell_data':dwell_dist_source,
        'cut_data':cut_dist_source,
        'pcut_data':post_dist_source,
        'unlooped_data':unlooped_dist_source,
        'loop_mat':ax_loop_mat,
        'dwell_mat':ax_dwell_mat,
        'percentile_source': perc_source,
        'percentile_view': perc_view,
        'cut_mat':ax_cut_mat}

# Define a reset button

reset_cb = CustomJS(args=args, code=reset_code + draw)
reset = Button(label="Click to reset plots, press ESC to clear selection")

reset.callback = reset_cb
loop_cb = CustomJS(args=args, code=loop_sel_code + sel_code + draw)
dwell_cb = CustomJS(args=args, code=dwell_sel_code + sel_code + draw)
cut_cb = CustomJS(args=args, code=cut_sel_code + sel_code + draw)

for a, t in zip([ax_loop_mat, ax_dwell_mat, ax_cut_mat], [loop_cb, dwell_cb, cut_cb]):
    tap_event = a.select(type=TapTool)
    tap_event.callback = t

endog_seq = vdj.io.endogenous_seqs()['WT12rss'][0]
for a, s  in zip([ax_loop_mat, ax_dwell_mat, ax_cut_mat], [loop_source, dwell_source, pcut_source]):
    a.yaxis.ticker = [0, 1, 2, 3]
    a.xaxis.ticker = np.arange(0, len(endog_seq), 1)
    ylab = {int(i):nt_idx[i] for i in range(4)}
    xlab = {int(i):b for i, b in zip(np.arange(0, len(endog_seq), 1), list(endog_seq))}
    a.yaxis.major_label_overrides = ylab
    a.xaxis.major_label_overrides = xlab

ax_loop.yaxis.visible = False


# Define titles
ax_loop_mat.title.text = "paired complex formation frequency"
ax_dwell_mat.title.text = "paired complex dwell time"
ax_cut_mat.title.text = "paired complex cleavage probability"

# Define the layout
spacer = Div(text="<br/>") #To give a little wiggle room between plots
leg_row = bokeh.layouts.row(ax_leg, reset)
loop_plots = bokeh.layouts.column(ax_loop_mat, ax_loop, spacer)
cut_plots = bokeh.layouts.column(ax_cut_mat, ax_cut, spacer)
dwell_row = bokeh.layouts.row(ax_dwell_unlooped, ax_dwell_cut, ax_dwell_all)    
col1 = bokeh.layouts.column(ax_loop_mat, ax_loop, ax_cut_mat, ax_cut,
ax_dwell_mat, dwell_row)

lay = bokeh.layouts.column(reset, ax_leg, ax_loop_mat, ax_leg2, ax_loop, ax_dwell_mat,
dwell_row, ax_cut_mat, ax_cut)

                           

# Define the color palettes
palette = bokeh.palettes.PRGn11
loop_color = LinearColorMapper(palette=palette, low=-0.25, high=0.25)
loop_bar = ColorBar(color_mapper=loop_color, location=(0, 0),
                    bar_line_color='black', ticker=FixedTicker(ticks=[-0.2, 0, 0.2]), 
                    width=15, height=50, background_fill_alpha=0)
ax_loop_mat.add_layout(loop_bar, 'right')

dwell_color = LinearColorMapper(palette=palette, low=-2.2, high=2.2)
dwell_bar = ColorBar(color_mapper=dwell_color, location=(0, 0),
                    bar_line_color='black', ticker=FixedTicker(ticks=[-2, 0, 2]), 
                    width=15, height=50, background_fill_alpha=0, title='[min]',
                    title_text_font_size='6pt')
ax_dwell_mat.add_layout(dwell_bar, 'right')

cut_color = LinearColorMapper(palette=palette, low=-0.62, high=0.62)
cut_bar = ColorBar(color_mapper=cut_color, location=(0, 0),
                    bar_line_color='black', ticker=FixedTicker(ticks=[-0.5, 0, 0.5]), 
                    width=15, height=50, background_fill_alpha=0) 
ax_cut_mat.add_layout(cut_bar, 'right')


# Define the color bars or the percentile
linear_mapper1 = LinearColorMapper(palette=colors1, low=5, high=99)
linear_mapper2 = LinearColorMapper(palette=colors2, low=5, high=99)
ticker = FixedTicker(ticks=[10, 25, 50, 75 ,95])
labels = {10:'10%', 25:'25%', 50:'50%', 75:'75%', 95:'95%'}
bar1 = ColorBar(color_mapper=linear_mapper1, ticker=ticker,
                location=(50, -10), border_line_color=None,
                major_label_overrides=labels, label_standoff=5, width=150, height=10,
                title='', background_fill_alpha=0, orientation='horizontal')
bar2 = ColorBar(color_mapper=linear_mapper2, ticker=ticker,
                location=(275, -10), border_line_color=None,
                major_label_overrides=labels, label_standoff=5, width=150, height=10,
                title='', background_fill_alpha=0, orientation='horizontal') 
ax_leg2.add_layout(bar1)
ax_leg2.add_layout(bar2)
ax_leg2.title.text = 'confidence interval'



# Populate the matrices.
loop_fig = ax_loop_mat.rect('pos', 'base_idx', width=1, height=1, source=loop_source, 
        fill_color=transform('diff', loop_color))
ax_loop_mat.tools.append(HoverTool(renderers=[loop_fig],
        tooltips=[('mutant', '@mutant'), 
                  ('difference in frequency', '@diff'),
                  ('looping frequency', '@loops_per_bead'),
                  ('number of beads', '@n_beads'),
                  ('number of loops', '@n_loops')]))
dwell_fig = ax_dwell_mat.rect('pos', 'base_idx', width=1, height=1, source=dwell_source, 
        fill_color=transform('diff', dwell_color))
ax_dwell_mat.tools.append(HoverTool(renderers=[dwell_fig],
        tooltips=[('mutant', '@mutant'), 
                  ('difference in median dwell time [min]', '@diff'),
                  ('median dwell time [min]', '@med_dwell'),
                  ('number of loops', '@n_loops')]))
cut_fig = ax_cut_mat.rect('pos', 'base_idx', width=1, height=1, source=pcut_source, 
        fill_color=transform('diff', cut_color))
ax_cut_mat.tools.append(HoverTool(renderers=[cut_fig],
        tooltips=[('mutant', '@mutant'), 
                  ('difference in cleavage probability', '@diff'),
                  ('cleavage probability', '@pcut'),
                  ('number of loops', '@n_loops'),
                  ('number of cuts', '@n_cuts')]))

# Add the reference features
ax_dwell_unlooped.step('dwell', 'ecdf',  color='slategrey', line_width=2,
alpha=1, source=unlooped_dist_ref)
# ax_dwell_unlooped.circle('dwell', 'ecdf',  size=4, fill_color='white', color='slategrey', alpha=1, source=unlooped_dist_ref)
ax_dwell_cut.step('dwell', 'ecdf',  color='slategrey', alpha=0.8,
source=cut_dist_ref, line_width=2)
# ax_dwell_cut.circle('dwell', 'ecdf',  size=4, fill_color='white', color='slategrey', alpha=1, source=cut_dist_ref)
ax_dwell_all.step('dwell', 'ecdf',  color='slategrey', alpha=0.8,
source=dwell_dist_ref, line_width=2)
# ax_dwell_all.circle('dwell', 'ecdf', size=4, fill_color='white', color='slategrey', alpha=1, source=dwell_dist_ref)


ax_cut.line('probability', 'posterior', source=post_ref, color='slategrey',
            alpha=1)
ax_cut.varea('probability', 0, 'posterior', source=post_ref, fill_color='slategrey',
            alpha=0.8)


# Add the point mutant features.
mut_rep = ax_loop.triangle('loops_per_bead', 0.5, color='dodgerblue', alpha=0.5, 
                 source=loop_source, view=loop_view, size=8)
 
ax_dwell_unlooped.step('dwell', 'ecdf',color='dodgerblue', alpha=1,
                 source=unlooped_dist_source, view=unlooped_view, line_width=2)
# ax_dwell_unlooped.circle('dwell', 'ecdf', size=4, color='dodgerblue',
                #  alpha=1, source=unlooped_dist_source, view=unlooped_view,
                #  fill_color='white')
ax_dwell_cut.step('dwell', 'ecdf',color='dodgerblue', alpha=1,
                 source=cut_dist_source, view=cut_view, line_width=2)
# ax_dwell_cut.circle('dwell', 'ecdf', size=4, fill_color='white', color='dodgerblue',
                #  alpha=0.8, source=cut_dist_source, view=cut_view)
ax_dwell_all.step('dwell', 'ecdf',color='dodgerblue', alpha=1,
                 source=dwell_dist_source, view=dwell_view, line_width=2)
# ax_dwell_all.circle('dwell', 'ecdf', size=4, fill_color='white', color='dodgerblue', alpha=0.8,
                #  source=dwell_dist_source, view=dwell_view)

ax_cut.line('probability', 'posterior', source=post_dist_source, color='dodgerblue',
            view=post_view)
ax_cut.varea('probability', y1=0, y2='posterior', source=post_dist_source, 
            fill_color='dodgerblue', view=post_view, alpha=0.8)



# Plot x's for virgin mutation
x, y = [], []
for p in range(28):
   for b in range(4):
      if (len(pcut_mat[(pcut_mat['pos'] == p) &
         (pcut_mat['base'] == nt_idx[b])]) == 0) & (ref_seq[p] != b):
          x.append(p)
          y.append(b)

for a in [ax_loop_mat, ax_dwell_mat, ax_cut_mat]:
    a.x(x, y, color='slategrey')
    # Fill in the wild-type positions
    a.rect(np.arange(0, 29, 1), ref_seq, width=1, height=1, fill_color='white',
            line_color='slategrey')
    a.circle(np.arange(0, 29, 1), ref_seq,  fill_color='#f5e3b3', alpha=0.75,
              line_color='slategrey', line_width=1, size=10)



# Adjust legend as necessary
ax_loop.legend.spacing = 1
ax_loop.legend.padding = 4


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
bokeh.plotting.output_file('./point_voyager.html', mode='inline')

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(lay)


#%%
