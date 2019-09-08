#%%
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, 
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool)
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
        'loops_per_bead':d['n_loops'].sum() / len(d['bead_idx']),
        'n_beads':len(d['bead_idx']), 'n_loops':d['n_loops'].sum()},
        ignore_index=True)
    for _g, _d in d.groupby(['replicate']):
        rep_df = rep_df.append({'mutant':g,
            'loops_per_bead':_d['n_loops'].sum() / len(_d['bead_idx']),
            'n_beads':len(_d['bead_idx']), 'n_loops':_d['n_loops'].sum()},
            ignore_index=True)
pooled_df['y'] = 0
rep_df['y'] = np.random.normal(0,0.05, len(rep_df))
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
for g, d in df.groupby(['mutant']):
    mut_seq = list(seq[g][0])
    for i, b in enumerate(mut_seq):
        if i not in invar_pos:
            var_dict['pos'].append(i)
            var_dict['base'].append(b)
            var_dict['y'].append(0)
            var_dict['mutant'].append(g)


#%%
# ##############################################################################
# SOURCE AND VIEW DEFINITION
# ##############################################################################
# Set up the dropdown with the mutations
menu = [ ("DFL 16.1-3'", 'DFL1613'), ("DFL 16.1-5'", 'DFL161'), 
         ('V1-135', 'V1-135'), ('V9-120', 'V9-120'), ('V10-96', 'V10-96'),
         ('V19-93', 'V19-93'), ('V4-57-1 (reference)', 'WT12rss'), 
         ('V4-55', 'V4-55'), ('V5-43', 'V5-43'), ('V8-18', 'V8-18'), 
         ('V6-17', 'V6-17'), ('V6-15', 'V6-15')]
mut_sel = Dropdown(value='V4-57-1 (reference)', menu=menu)

# Define the filter on the mutant props.
mut_filter = GroupFilter(column_name="mutant", group=mut_sel.value)

# Define the sources
variant = ColumnDataSource(var_dict)
dwell_source = ColumnDataSource(dwell_hist)
cut_source = ColumnDataSource(cut_hist)
post_source = ColumnDataSource(posteriors)
pooled_source = ColumnDataSource(pooled_df)
rep_source = ColumnDataSource(rep_df)

# Define the Views
seq_view = CDSView(source=variant, filters=[mut_filter])
dwell_view = CDSView(source=dwell_source, filters=[mut_filter])
cut_view = CDSView(source=cut_source, filters=[mut_filter])
pooled_loop_view = CDSView(source=pooled_source, filters=[mut_filter])
rep_loop_view = CDSView(source=rep_source, filters=[mut_filter])
post_view = CDSView(source=post_source, filters=[mut_filter])


#%% 
# ##############################################################################
# DEFINE THE CANVASES
# ##############################################################################
ax_seq = bokeh.plotting.figure(height=40, y_range=[0, 0.1], tools=[''])   
ax_loop = bokeh.plotting.figure(height=120, 
                                x_axis_label='paired complexes per bead',
                                x_range=[0, 0.8], y_range=[-0.5, 0.5], tools=[''])
ax_dwell = bokeh.plotting.figure(height=200, x_axis_label='paired complex dwell time [min]',
                                y_axis_label='number of observations', tools=[''])
ax_cut = bokeh.plotting.figure(height=200, x_axis_label='cutting probability',
                               y_axis_label='posterior probability', tools=[''])

# Set features of the plots
ax_seq.xaxis.visible = False
ax_seq.yaxis.visible = False
ax_seq.grid.visible = False
ax_loop.yaxis.visible = False

# Add hover tooltips to the loop plot 
tooltips = [('# beads', '@n_beads'), ('# paired complexes', '@n_loops')]
ax_loop.add_tools(HoverTool(tooltips=tooltips))

# Define the layout
lay = bokeh.layouts.column(mut_sel, ax_seq, ax_loop, ax_dwell, ax_cut) 

# ##############################################################################
# POPULATE CANVASES
# ############################################################################## 
# Sequence
invar_glyph = bokeh.models.glyphs.Text(x='pos', y='y', text='base', text_font='Courier',
                                    text_color='#95a3b2', text_font_size='28px')
var_glyph = bokeh.models.glyphs.Text(x='pos', y='y', text='base', text_font='Courier',
                                    text_color='dodgerblue', text_font_size='28px', 
                                    text_alpha=0.8)
ax_seq.add_glyph(invariant, invar_glyph)
ax_seq.add_glyph(variant, var_glyph, view=seq_view)

# Loops per bead
ax_loop.x(x='loops_per_bead', y='y', source=rep_source,
                view=rep_loop_view, color='slategray', legend='replicate',
                alpha=0.5, size=8)
ax_loop.circle(x='loops_per_bead', y='y', source=pooled_source,
                view=pooled_loop_view, color='dodgerblue', legend='pooled',
                size=10, fill_alpha=0.5)

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
# Set up the callback
args = {'sel':mut_sel, 'filter':mut_filter,
        'seq_view':seq_view, 'pooled_view':pooled_loop_view, 
        'rep_view':rep_loop_view, 'dwell_view':dwell_view, 'cut_view':cut_view,
        'post_view':post_view,
        'seq_data':variant, 'pooled_data':pooled_source, 'rep_data':rep_source,
        'dwell_data':dwell_source, 'cut_data':cut_source, 'post_data':post_source}
callback = CustomJS(args=args, code="""
                    var mut = sel.value;
                    filter.group = mut;
                    var views = [seq_view, pooled_view, rep_view, dwell_view, 
                                cut_view, post_view];
                    var data = [seq_data, pooled_data, rep_data, dwell_data, 
                                cut_data, post_data];
                    for (var i = 0; i < views.length; i++) { 
                          views[i].filters[0] = filter;
                          data[i].data.view = views[i];
                         data[i].change.emit();}
                    """)

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



#%%
