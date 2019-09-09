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
var_dict = {'pos':[], 'base':[], 'mutant':[], 'y1':[], 'y2':[]}
for g, d in df.groupby(['mutant']):
    mut_seq = list(seq[g][0])
    for i, b in enumerate(mut_seq):
        if i not in invar_pos:
            var_dict['pos'].append(i)
            var_dict['base'].append(b)
            var_dict['y1'].append(1)
            var_dict['y2'].append(2)
            var_dict['mutant'].append(g)



# ##############################################################################
# SOURCE AND VIEW DEFINITION
# ##############################################################################
# Set up the dropdown with the mutations
menu = [ ("DFL16.1-3'", 'DFL1613'), ("DFL16.1-5'", 'DFL161'), 
         ('V1-135', 'V1-135'), ('V9-120', 'V9-120'), ('V10-96', 'V10-96'),
         ('V19-93', 'V19-93'), ('V4-57-1 (reference)', 'WT12rss'), 
         ('V4-55', 'V4-55'), ('V5-43', 'V5-43'), ('V8-18', 'V8-18'), 
         ('V6-17', 'V6-17'), ('V6-15', 'V6-15'), ('None ', 'nan')]
mut_sel1 = Dropdown(label='Select Endogenous 12RSS', value="None", 
                    menu=menu, button_type='primary')
mut_sel2 = Dropdown(label='Select Endogenous 12RSS', value="None", 
                    menu=menu, button_type='danger')

# Define the filter on the mutant props.
mut_filter1 = GroupFilter(column_name="mutant", group=mut_sel1.value)
mut_filter2 = GroupFilter(column_name="mutant", group=mut_sel2.value)

# Define the sources
variant1 = ColumnDataSource(var_dict)
dwell_source1 = ColumnDataSource(dwell_hist)
cut_source1 = ColumnDataSource(cut_hist)
post_source1 = ColumnDataSource(posteriors)
pooled_source1 = ColumnDataSource(pooled_df)
rep_source1 = ColumnDataSource(rep_df)
variant2 = ColumnDataSource(var_dict)
dwell_source2 = ColumnDataSource(dwell_hist)
cut_source2 = ColumnDataSource(cut_hist)
post_source2 = ColumnDataSource(posteriors)
pooled_source2 = ColumnDataSource(pooled_df)
rep_source2 = ColumnDataSource(rep_df)

# Define the Views
seq_view1 = CDSView(source=variant1, filters=[mut_filter1])
seq_view2 = CDSView(source=variant2, filters=[mut_filter2])
dwell_view1 = CDSView(source=dwell_source1, filters=[mut_filter1])
dwell_view2 = CDSView(source=dwell_source2, filters=[mut_filter2])
cut_view1 = CDSView(source=cut_source1, filters=[mut_filter1])
cut_view2 = CDSView(source=cut_source2, filters=[mut_filter2])
pooled_loop_view1 = CDSView(source=pooled_source1, filters=[mut_filter1])
pooled_loop_view2 = CDSView(source=pooled_source2, filters=[mut_filter2])
rep_loop_view1 = CDSView(source=rep_source1, filters=[mut_filter1])
rep_loop_view2 = CDSView(source=rep_source2, filters=[mut_filter2])
post_view1 = CDSView(source=post_source1, filters=[mut_filter1])
post_view2 = CDSView(source=post_source2, filters=[mut_filter2])


#%% 
# ##############################################################################
# DEFINE THE CANVASES
# ##############################################################################
ax_seq = bokeh.plotting.figure(height=100, y_range=[0, 4], tools=[''])   
ax_loop = bokeh.plotting.figure(height=140, 
                                x_axis_label='paired complexes per bead\n   ',
                                x_range=[-0.05, 0.95], y_range=[-0.5, 0.5], tools=[''])
ax_dwell = bokeh.plotting.figure(height=200, x_axis_label='paired complex dwell time [min]\n  ',
                                y_axis_label='number of observations', tools=[''])
ax_cut = bokeh.plotting.figure(height=200, x_axis_label='cutting probability',
                               y_axis_label='posterior probability', tools=[''])

description1 = Div(text="""<center>Select an RSS</center>""")
description2 = Div(text="""<center>Select an RSS</center>""")

for a in [ax_seq, ax_loop, ax_dwell, ax_cut]:
    a.toolbar.logo = None
# Set features of the plots
ax_seq.xaxis.visible = False
ax_seq.yaxis.visible = False 
ax_seq.xgrid.visible = False
ax_seq.ygrid.visible = False
ax_seq.background_fill_color = None
ax_seq.border_fill_color = None
ax_seq.outline_line_color = None
ax_seq.grid.visible = False
ax_seq.title.text_font_style = 'bold'
ax_loop.yaxis.visible = False

# Add titles
ax_loop.title.text = 'paired complex formation frequency'
ax_dwell.title.text = 'paired complex dwell time distribution'
ax_cut.title.text = 'paired complex cleavage probability'

# Add hover tooltips to the loop plot 
tooltips = [('mutant', '@mutant'), ('# beads', '@n_beads'), ('# paired complexes', '@n_loops')]
ax_loop.add_tools(HoverTool(tooltips=tooltips))

# Define the layout
selections = bokeh.layouts.row(mut_sel1, mut_sel2)
descriptions = bokeh.layouts.row(description1, description2)
lay = bokeh.layouts.column(selections, descriptions, bokeh.layouts.Spacer(), ax_loop, ax_dwell, ax_cut)
lay.sizing_mode = 'scale_width'
# ############################################################################## # POPULATE CANVASES
# ############################################################################## 
# Sequence
invar_glyph = bokeh.models.glyphs.Text(x='pos', y='y', text='base', text_font='Courier',
                                    text_color='#95a3b2', text_font_size='28px')
var_glyph1 = bokeh.models.glyphs.Text(x='pos', y='y1', text='base', text_font='Courier',
                                    text_color='dodgerblue', text_font_size='28px', 
                                    text_alpha=0.8)
var_glyph2 = bokeh.models.glyphs.Text(x='pos', y='y2', text='base', text_font='Courier',
                                    text_color='tomato', text_font_size='28px', 
                                    text_alpha=0.8)
ax_seq.add_glyph(invariant, invar_glyph)
ax_seq.add_glyph(variant1, var_glyph1, view=seq_view1)
ax_seq.add_glyph(variant2, var_glyph2, view=seq_view2)

# Loops per bead
ax_loop.triangle(x='loops_per_bead', y='y', source=rep_source1,
                view=rep_loop_view1, color='dodgerblue', legend='replicate',
                alpha=0.5, size=8)
ax_loop.triangle(x='loops_per_bead', y='y', source=rep_source2,
                view=rep_loop_view2, color='tomato', legend='replicate',
                alpha=0.5, size=8)

ax_loop.circle(x='loops_per_bead', y='y', source=pooled_source1,
                view=pooled_loop_view1, color='dodgerblue', legend='pooled',
                size=10, fill_alpha=0.5)

ax_loop.circle(x='loops_per_bead', y='y', source=pooled_source2,
                view=pooled_loop_view2, color='tomato', legend='pooled',
                size=10, fill_alpha=0.5)
ax_loop.legend.click_policy = 'hide'

# Add the histogram
ax_dwell.quad(left='left', bottom='bottom', top='top', right='right', 
               view=dwell_view1, source=dwell_source1, color='dodgerblue',
               alpha=0.5, legend="all PC events")
ax_dwell.quad(left='left', bottom='bottom', top='top', right='right', 
               view=dwell_view2, source=dwell_source2, line_color='tomato',
               alpha=0.5, fill_color=None, line_width=2.5, legend="all PC events")
ax_dwell.quad(left='left', bottom='bottom', top='top', right='right', 
               view=cut_view1, source=cut_source1, hatch_color='navy',
               alpha=0.5, hatch_pattern='/', legend="cleavage events")
ax_dwell.quad(left='left', bottom='bottom', top='top', right='right', 
               view=cut_view2, source=cut_source2, fill_color=None, hatch_color='firebrick',
               alpha=0.5, hatch_pattern='\\', legend="cleavage events")
ax_dwell.legend.click_policy = 'hide'
# Cutting probability posterior
ax_cut.line(x='probability', y='posterior', source=post_source1, view=post_view1, 
            color='dodgerblue')
ax_cut.varea(x='probability', y1=0, y2='posterior', fill_color='dodgerblue',
            fill_alpha=0.5, source=post_source1, view=post_view1)
ax_cut.line(x='probability', y='posterior', source=post_source2, view=post_view2, 
            color='tomato')
ax_cut.varea(x='probability', y1=0, y2='posterior', fill_color='tomato',
            fill_alpha=0.5, source=post_source2, view=post_view2)

# ##############################################################################
# CALLBACK DEFINITION
# ##############################################################################
# Set up the callback
args = {'sel1':mut_sel1, 'filter1':mut_filter1, 'sel2':mut_sel2, 'filter2':mut_filter2,
        'seq_view1':seq_view1, 'seq_view2':seq_view2, 
        'pooled_view1':pooled_loop_view1, 'pooled_view2':pooled_loop_view2, 
        'rep_view1':rep_loop_view1, 'rep_view2':rep_loop_view2,
        'dwell_view1':dwell_view1, 'dwell_view2':dwell_view2,
        'cut_view1':cut_view1, 'cut_view2':cut_view2,
        'post_view1':post_view1, 'post_view2':post_view2,
        'seq_data1':variant1,  'seq_data2':variant2, 
        'pooled_data1':pooled_source1, 'pooled_data2':pooled_source2, 
        'rep_data1':rep_source1, 'rep_data2':rep_source2,
        'dwell_data1':dwell_source1, 'dwell_data2':dwell_source2, 
        'cut_data1':cut_source1, 'cut_data2':cut_source2, 
        'post_data1':post_source1, 'post_data2':post_source2,
        'description1':description1, 'description2':description2}

callback = CustomJS(args=args, code="""
  var mut1 = sel1.value;
  var mut2 = sel2.value;
  filter1.group = mut1;
  filter2.group = mut2;
  var views1 = [seq_view1, pooled_view1, rep_view1, dwell_view1, cut_view1, post_view1]; 
  var views2 = [seq_view2, pooled_view2, rep_view2, dwell_view2, cut_view2, post_view2];
  var data1 = [seq_data1, pooled_data1, rep_data1, dwell_data1, cut_data1, post_data1];
  var data2 = [seq_data2, pooled_data2, rep_data2, dwell_data2, cut_data2, post_data2];
  for (var i = 0; i < views1.length; i++) { 
       views1[i].filters[0] = filter1;
       data1[i].data.view = views1[i];
       data1[i].change.emit();}
  for (var i = 0; i < views2.length; i++) { 
       views2[i].filters[0] = filter2;
       data2[i].data.view = views2[i];
       data2[i].change.emit();}

 // Set the descriptions of the sequences 
if (mut1 === 'WT12rss') {
     var desc1 = 'V4-57-1 (reference)';
  }
else if (mut1 === 'DFL161') {
    var desc1 = "DFL16.1-5'";
}
else if (mut1 === 'DFL1613') {
    var desc1 = "DFL16.1-3'";
}
else if (mut1 === 'nan') { 
    var desc1 = "None Selected"; 
}
else { var desc1 = mut1}

if (mut2 === 'WT12rss') {
     var desc2 = 'V4-57-1 (reference)';
  }
else if (mut2 === 'DFL161') {
    var desc2 = "DFL16.1-5'";
}
else if (mut2 === 'DFL1613') {
    var desc1 = "DFL16.1-3'";
}
else if (mut2 === 'nan') {
    var  desc2 = "None Selected";
}
else { var desc2 = mut2}

  description1.text = "<span style='color: dodgerblue; font-size: 14pt; width: 50% margin: 0 auto;'>      " + desc1 + "    </span>"
  description2.text = "<span style='color: tomato; font-size: 14pt;'>      " + desc2 + "     </span>"
        """)

mut_sel1.js_on_change('value', callback)
mut_sel2.js_on_change('value', callback)


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
                'text_font_style': 'bold',
                'align': 'center',
                'text_font': 'Helvetica',

                'offset': 2,
            }}}

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(lay)



#%%
