#%% 
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, 
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool,
                          Segment, ColorBar, LinearColorMapper, FixedTicker)
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.embed import components
import bokeh.palettes
import vdj.io
import vdj.stats
bokeh.plotting.output_file('../../figures/interactives/endogenous_sequence_exporer.html',
                        mode='inline')


#%%
# ##############################################################################
# LOAD DATA SETS AND RESTRICT TO ENDOGENOUS
# ##############################################################################
# Load the necessary data sets
dwell_times = pd.read_csv('../../data/compiled_dwell_times.csv', comment='#')
posteriors = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv', comment='#')
loops = pd.read_csv('../../data/compiled_looping_frequency_bootstrap.csv', comment='#')

# Identify the endogenous muts
dfs  = []
for i, df in enumerate([dwell_times, posteriors, loops]):
    # Rename the point mutant 12SpacC1A to  the endogenous mutant V4-55
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
unloop_dwells = dwell_times[dwell_times['cut']==0]




#%%
# ##############################################################################
# GENERATE HISTOGRAMS OF DWELL TIMES 
# ##############################################################################
# Generate the histogrammed dwell times.
dfs = []
for source in [dwell_times, cut_dwells, unloop_dwells]: 
    bin_dfs = []
    for g, d in source.groupby('mutant'):
        x, y = np.sort(d['dwell_time_min'].values), np.arange(0, len(d), 1) / len(d)
        y[-1] = 1
        # hist, bins = np.histogram(d['dwell_time_min'], bins=bins)
        _df = pd.DataFrame()
        _df['dwell_time'] = x
        _df['ecdf'] = y
        _df['mutant'] = g
        bin_dfs.append(_df)
    dwell_hist = pd.concat(bin_dfs)
    dfs.append(dwell_hist)
dwell_dist, cut_dist, unloop_dist = dfs


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
                    menu=menu)
mut_sel2 = Dropdown(label='Select Endogenous 12RSS', value="None", 
                    menu=menu, button_type='primary')

# Define the filter on the mutant props.
mut_filter1 = GroupFilter(column_name="mutant", group=mut_sel1.value)
mut_filter2 = GroupFilter(column_name="mutant", group=mut_sel2.value)

# Define the sources
leg_source = ColumnDataSource({'xs':[[-10, -1], [-10, -1]], 
                               'ys':[[-10, -1], [-10, -1]], 
                               'c':['slategrey', 'dodgerblue'], 'mutant':['None Selected', 'None Selected']})

dwell_all_source1 = ColumnDataSource(dwell_dist)
dwell_cut_source1 = ColumnDataSource(cut_dist)
dwell_unloop_source1 = ColumnDataSource(unloop_dist)
post_source1 = ColumnDataSource(posteriors)
loop_source1 = ColumnDataSource(loops)
dwell_all_source2 = ColumnDataSource(dwell_dist)
dwell_cut_source2 = ColumnDataSource(cut_dist)
dwell_unloop_source2 = ColumnDataSource(unloop_dist)
cut_source2 = ColumnDataSource(cut_dist)
post_source2 = ColumnDataSource(posteriors)
loop_source2 = ColumnDataSource(loops)

# Define the Views
dwell_all_view1 = CDSView(source=dwell_all_source1, filters=[mut_filter1])
dwell_all_view2 = CDSView(source=dwell_all_source2, filters=[mut_filter2])
dwell_cut_view1 = CDSView(source=dwell_cut_source1, filters=[mut_filter1])
dwell_cut_view2 = CDSView(source=dwell_cut_source2, filters=[mut_filter2])
dwell_unloop_view1 = CDSView(source=dwell_unloop_source1, filters=[mut_filter1])
dwell_unloop_view2 = CDSView(source=dwell_unloop_source2, filters=[mut_filter2])
post_view1 = CDSView(source=post_source1, filters=[mut_filter1])
post_view2 = CDSView(source=post_source2, filters=[mut_filter2])
loop_view1 = CDSView(source=loop_source1, filters=[mut_filter1])
loop_view2 = CDSView(source=loop_source2, filters=[mut_filter2])


#%% 
# ##############################################################################
# DEFINE THE CANVASES
# ##############################################################################
ax_leg = bokeh.plotting.figure(height=80, width=600, tools =[''], toolbar_location=None, x_range=[0, 1], y_range=[0, 1])

ax_loop = bokeh.plotting.figure(height=140, width=800, 
                                x_axis_label='paired complexes per bead\n   ',
                                x_range=[-0.05, 0.95], y_range=[-0.8, 0.8], tools=[''])
ax_dwell_all = bokeh.plotting.figure(height=300, width=400, x_axis_label='paired complex dwell time [min]',
                                y_axis_label='ECDF', tools=[''], x_axis_type='log')
ax_dwell_cut = bokeh.plotting.figure(height=300, width=400,x_axis_label='paired complex dwell time [min]',
                                y_axis_label='ECDF', tools=[''], x_axis_type='log')
ax_dwell_unloop = bokeh.plotting.figure(height=300, width=400, x_axis_label='paired complex dwell time [min]',
                                y_axis_label='ECDF', tools=[''], x_axis_type='log')

ax_cut = bokeh.plotting.figure(height=300, width=400, x_axis_label='cutting probability',
                               y_axis_label='posterior probability', tools=[''])

# Add legend entries
ax_leg.multi_line('xs', 'ys', line_width=10, color='c', legend='mutant', source=leg_source)
ax_leg.legend.location = 'center'
ax_leg.legend.orientation  = 'horizontal'
ax_leg.background_fill_color = None
ax_leg.outline_line_color = None
ax_leg.legend.background_fill_color = None
ax_leg.xaxis.visible = False
ax_leg.yaxis.visible = False

for a in [ax_loop, ax_dwell_all, ax_dwell_cut, ax_dwell_unloop, ax_cut]:
    a.toolbar.logo = None
# Set features of the plots
ax_loop.yaxis.visible = False

# Add titles
ax_loop.title.text = 'paired complex formation frequency'
ax_dwell_unloop.title.text = 'unlooped PCs'
ax_dwell_cut.title.text = 'cleaved PCs'
ax_dwell_all.title.text = 'all PCs'
ax_cut.title.text = 'paired complex cleavage probability'

# Define the layout
selections = bokeh.layouts.row(mut_sel1, mut_sel2)
dwell_row = bokeh.layouts.row(ax_dwell_unloop, ax_dwell_cut)
final_row = bokeh.layouts.row(ax_dwell_all, ax_cut)
col = bokeh.layouts.column(ax_loop, dwell_row, final_row)
spacer = Div(text="<br/>")
lay = bokeh.layouts.column(selections, ax_leg, col)
# lay.sizing_mode = 'scale_width'
# ############################################################################## # POPULATE CANVASES
# ############################################################################## 
# Loops per bead
colors1 = bokeh.palettes.Greys9[1:-2]
colors2 = bokeh.palettes.Blues9[1:-2]


# Define colorbars
linear_mapper1 = LinearColorMapper(palette=colors1, low=5, high=99)
linear_mapper2 = LinearColorMapper(palette=colors2, low=5, high=99)
ticker = FixedTicker(ticks=[10, 25, 50, 75 ,95])
labels = {10:'10%', 25:'25%', 50:'50%', 75:'75%', 95:'95%'}
bar1 = ColorBar(color_mapper=linear_mapper1, ticker=ticker,
                location=(0, 0), border_line_color=None,
                major_label_overrides=labels, label_standoff=10, width=150, orientation='horizontal',
                title='confidence interval')
bar2 = ColorBar(color_mapper=linear_mapper2, ticker=ticker,
                location=(0, 0), border_line_color=None,
                major_label_overrides=labels, label_standoff=10, width=150, orientation='horizontal',
                title='confidence interval') 
ax_leg.add_layout(bar1, 'left')
ax_leg.add_layout(bar2, 'right')


percentile_source1 = []
perc_view1 = []
perc_view2 = []
percentile_source2 = []

percs = list(np.sort(loops['percentile'].unique()))
percs.reverse()

mut1_circ = ax_loop.triangle(name='mut1', x='loops_per_bead', y=0.5, source=loop_source1,
                view=loop_view1, fill_color='white', line_color='slategrey', 
                size=10, level='overlay', legend='observed frequency')
mut2_circ = ax_loop.triangle(name='mut2', x='loops_per_bead', y=-0.5, source=loop_source2,
                view=loop_view2, fill_color='white', line_color='dodgerblue', 
                size=10, level='overlay', legend='observed frequency')

for i, p in enumerate(percs):
    d = loops[loops['percentile']==p]
    _source1 = ColumnDataSource(d)
    _source2 = ColumnDataSource(d)
    _view1 = CDSView(source=_source1, filters=[mut_filter1])
    _view2 = CDSView(source=_source2, filters=[mut_filter2])
    percentile_source1.append(_source1)
    percentile_source2.append(_source2)
    perc_view1.append(_view1)
    perc_view2.append(_view2)

    band1 = Segment(x0='low', x1='high', y0=0.5, y1=0.5, 
                    line_color=colors1[-1 * (i+1)], line_width=25) 
    band2 = Segment(x0='low', x1='high', y0=-0.5, y1=-0.5, 
                    line_color=colors2[-1 * (i+1)], line_width=25) 
        
    ax_loop.add_glyph(_source1, band1, view=_view1)
    ax_loop.add_glyph(_source2, band2, view=_view2)


# Add the histogram
ax_dwell_all.step('dwell_time', 'ecdf', source=dwell_all_source1, view=dwell_all_view1,
                    color='slategrey', line_width=2)
ax_dwell_all.step('dwell_time', 'ecdf', source=dwell_all_source2, view=dwell_all_view2,
                    color='dodgerblue', line_width=2)
ax_dwell_cut.step('dwell_time', 'ecdf', source=dwell_cut_source1, view=dwell_cut_view1,
                    color='slategrey', line_width=2)
ax_dwell_cut.step('dwell_time', 'ecdf', source=dwell_cut_source2, view=dwell_cut_view2,
                    color='dodgerblue', line_width=2)
ax_dwell_unloop.step('dwell_time', 'ecdf', source=dwell_unloop_source1, view=dwell_unloop_view1,
                    color='slategrey', line_width=2)
ax_dwell_unloop.step('dwell_time', 'ecdf', source=dwell_unloop_source2, view=dwell_unloop_view2,
                    color='dodgerblue', line_width=2)




# Cutting probability posterior
ax_cut.line(x='probability', y='posterior', source=post_source1, view=post_view1, 
            color='slategrey')
ax_cut.varea(x='probability', y1=0, y2='posterior', fill_color='slategrey',
            fill_alpha=0.5, source=post_source1, view=post_view1)
ax_cut.line(x='probability', y='posterior', source=post_source2, view=post_view2, 
            color='dodgerblue')
ax_cut.varea(x='probability', y1=0, y2='posterior', fill_color='dodgerblue',
            fill_alpha=0.5, source=post_source2, view=post_view2)

# ##############################################################################
# CALLBACK DEFINITION
# ##############################################################################
# Set up the callback
args = {'sel1':mut_sel1, 'filter1':mut_filter1, 'sel2':mut_sel2, 'filter2':mut_filter2,
        'loop_view1':loop_view1, 'loop_view2':loop_view2, 
        'dwell_all_view1':dwell_all_view1, 'dwell_all_view2':dwell_all_view2,
        'dwell_cut_view1':dwell_cut_view1, 'dwell_cut_view2':dwell_cut_view2,
        'dwell_unloop_view1':dwell_unloop_view1, 'dwell_unloop_view2':dwell_unloop_view2,
        'post_view1':post_view1, 'post_view2':post_view2,
        'loop_data1':loop_source1, 'loop_data2':loop_source2,
        'dwell_cut_data1':dwell_cut_source1, 'dwell_cut_data2':dwell_cut_source2, 
        'dwell_unloop_data1':dwell_unloop_source1, 'dwell_unloop_data2':dwell_unloop_source2, 
        'dwell_all_data1':dwell_all_source1, 'dwell_all_data2':dwell_all_source2, 
        'post_data1':post_source1, 'post_data2':post_source2,
        'percentile_source1': percentile_source1,
        'percentile_source2': percentile_source2,
        'percentile_view1': perc_view1,
        'percentile_view2': perc_view2,
        'leg_source':leg_source}

callback = CustomJS(args=args, code="""
  var mut1 = sel1.value;
  var mut2 = sel2.value;
  filter1.group = mut1;
  filter2.group = mut2;

  // Iterate through each percentile. 
  for (var i = 0; i < percentile_source1.length; i++ ) {
      var perc1 = percentile_source1[i];
      var perc2 = percentile_source2[i];
      percentile_view1[i].filters[0] = filter1;
      percentile_view2[i].filters[0] = filter2;
      perc1.data.view = percentile_view1[i];
      perc2.data.view = percentile_view2[i];
      perc1.change.emit();
      perc2.change.emit();

  }

  var views1 = [loop_view1, dwell_all_view1, dwell_cut_view1, dwell_unloop_view1, post_view1]; 
  var views2 = [loop_view2, dwell_all_view2, dwell_cut_view2, dwell_unloop_view2, post_view2];
  var data1 = [loop_data1, dwell_all_data1, dwell_cut_data1, dwell_unloop_data1, post_data1];
  var data2 = [loop_data2, dwell_all_data2, dwell_cut_data2, dwell_unloop_data2, post_data2];
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
    var desc2 = "DFL16.1-3'";
}
else if (mut2 === 'nan') {
    var  desc2 = "None Selected";
}
else { var desc2 = mut2}
    leg_source.data['mutant'] = [desc1, desc2];
    leg_source.change.emit()
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
