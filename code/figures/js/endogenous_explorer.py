#%%
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, 
                          CDSView, Grid, GroupFilter, Band, Dropdown)
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.embed import components
import vdj.io
import vdj.stats
bokeh.plotting.output_file('./endogenous_voyager.html')


#%%
# ##############################################################################
# LOAD DATA SETS AND RESTRICT TO ENDOGENOUS
# ##############################################################################
# Load the necessary data sets
dwell_times = pd.read_csv('../../../data/compiled_dwell_times.csv')
posteriors = pd.read_csv('../../../data/pooled_cutting_probability_posteriors.csv')
loops = pd.read_csv('../../../data/compiled_looping_events.csv')

# Identify the endogenous muts
dfs  = []
for i, df in enumerate([dwell_times, posteriors, loops]):
    for g, d in df.groupby(['mutant']):
        if ('12Spac' not in g) & ('12Hept' not in g) & ('12Non' not in g):
            mut_class = 'endogenous'
        else:
            mut_class = 'point' 
       
        df.loc[df['mutant']==g, 'class'] = mut_class
    df = df[(df['class']=='endogenous') & (df['salt']=='Mg') & (df['hmgb1']==80)]
    dfs.append(df)
dwell_times, posteriors, loops = dfs
cut_dwells = dwell_times[dwell_times['cut']==1]


# %%
# ##############################################################################
# COMPUTE POOLED v REPLICATE LOOP FREQUENCIES
# ##############################################################################
pooled_df = pd.DataFrame()
rep_df = pd.DataFrame()
for g, d in loops.groupby('mutant'):
    pooled_df = pooled_df.append({'mutant': g, 
        'loops_per_bead':d['n_loops'].sum() / len(d['bead_idx'])},
        ignore_index=True)
    for _g, _d in d.groupby(['replicate']):
        rep_df = rep_df.append({'mutant':g, 
            'loops_per_bead':d['n_loops'].sum() / len(d['bead_idx'])},
            ignore_index=True)
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

# %%
# ##############################################################################
#  DEFINE INVARIANT / VARIANT SEQUENCE MARKER
# ##############################################################################
# Find the invariable positions in the sequence
ind = np.zeros((28, 5))
seq = vdj.io.endogenous_seqs()
nt_idx = vdj.io.nucleotide_idx()

for g, d in df.groupby('mutant'):
    mut_seq = list(seq[g][0])
    for i, b, in enumerate(mut_seq):
        ind[i][nt_idx[b]] += 1

invar_pos = np.where(ind == len(df['mutant'].unique()))[0]
invar_base = [nt_idx[i] for i in np.where(ind == len(df['mutant'].unique()))[1]]
invariant = ColumnDataSource(dict(pos=invar_pos, 
                                  y= np.zeros(len(df['mutant'].unique())),
                                  base=invar_base))

# Set up a source for the variant positions
var_dict = {'pos':[], 'base':[], 'mutant':[], 'y':[]}
seq_view = ColumnDataSource(var_dict)
for g, d in df.groupby(['mutant']):
    mut_seq = list(seq[g][0])
    for i, b in enumerate(mut_seq):
        if i not in invar_pos:
            var_dict['pos'].append(i)
            var_dict['base'].append(b)
            var_dict['y'].append(0)
            var_dict['mutant'].append(g)

variant = ColumnDataSource(var_dict)

#%%
# ##############################################################################
# SOURCE AND VIEW DEFINITION
# ##############################################################################
# Set up the dropdown with the mutations
mut_sel = Select(value='', options=list(dwell_times['mutant'].unique()))

# Define the filter on the mutant props.
mut_filter = GroupFilter(column_name="mutant", group=mut_sel.value)

# Define the sources
dwell_source = ColumnDataSource(dwell_hist)
cut_source = ColumnDataSource(cut_hist)
post_source = ColumnDataSource(posteriors)
pooled_loop_source = ColumnDataSource(pooled_df)
rep_loop_source = ColumnDataSource(rep_df)

# Define the Views
dwell_view = CDSView(source=dwell_source, filters=[mut_filter])
cut_view = CDSView(source=cut_source, filters=[mut_filter])
pooled_loop_view = CDSView(source=pooled_loop_source, filters=[mut_filter])
rep_loop_view = CDSView(source=rep_loop_source, filters=[mut_filter])
seq_view = CDSView(source=variant, filters=[mut_filter])
post_view = CDSView(source=post_source, filters=[mut_filter])

#%% 
# ##############################################################################
# DEFINE THE CANVASES
# ##############################################################################
ax_seq = bokeh.plotting.figure(height=40, y_range=[0, 0.1])
ax_loop = bokeh.plotting.figure(height=100, 
                                x_axis_label='paired complexes per bead')
ax_dwell = bokeh.plotting.figure(height=200, x_axis_label='paired complex dwell time [min]',
                                y_axis_label='number of observations')
ax_cut = bokeh.plotting.figure(height=200, x_axis_label='cutting probability',
                               y_axis_label='posterior probability')

# Set features of the plots
ax_seq.xaxis.visible = False
ax_seq.yaxis.visible = False
ax_seq.grid.visible = False

# Define the layout
lay = bokeh.layouts.column(mut_sel, ax_seq, ax_loop, ax_dwell, ax_cut) 

# ##############################################################################
# POPULATE CANVASES
# ############################################################################## 
# Sequence
invar_glyph = bokeh.models.glyphs.Text(x='pos', y='y', text='base', text_font='Courier',
                                    text_color='#95a3b2', text_font_size='28px')
var_glyph = bokeh.models.glyphs.Text(x='pos', y='y', text='base', text_font='Courier',
                                    text_color='skyblue', text_font_size='28px')
ax_seq.add_glyph(invariant, invar_glyph)
ax_seq.add_glyph(variant, var_glyph, view=seq_view)

# Loops per bead
ax_loop.circle(x='loops_per_bead', y=0, source=rep_loop_source,
                view=rep_loop_view, color='slategray', legend='replicate')
ax_loop.circle(x='loops_per_bead', y=0, source=pooled_loop_source,
                view=pooled_loop_view, color='dodgerblue', legend='pooled')

# Add the histogram
ax_dwell.quad(left='left', bottom='bottom', top='top', right='right', 
               view=dwell_view, source=dwell_source, color='dodgerblue',
               alpha=0.5, legend="all PC events")
ax_dwell.quad(left='left', bottom='bottom', top='top', right='right', 
               view=cut_view, source=cut_source, color='slategrey',
               alpha=0.5, hatch_pattern='/', legend="cleavage events")

# Cutting probability posterior
ax_cut.line(x='probability', y='posterior', source=post_source, view=post_view, 
            color='dodgerblue')
ax_cut.varea(x='probability', y1=0, y2='posterior', fill_color='dodgerblue',
            fill_alpha=0.5, source=post_source, view=post_view)
#
# ##############################################################################
# CALLBACK DEFINITION
# ##############################################################################
with open('endogenous_explorer.js') as f:
    cb = f.read()

# Set up the callback
callback = CustomJS(code=cb, args={'sel':mut_sel, 'seq_view':seq_view, 'seqs':variant,
                          'filter':mut_filter, 'dwell_view':dwell_view, 'cut_view':cut_view,
                          'cut_hist':cut_source,
                          'dwell_hist':dwell_source, 'post_view':post_view,
                          'post':post_source})
mut_sel.js_on_change('value', callback)


#%% 
# ##############################################################################
# THEMING AND FILE I/O
# ##############################################################################
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
bokeh.io.save(lay)

