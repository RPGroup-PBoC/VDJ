#%% 
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, Text,
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool,
                          IndexFilter, TapTool)
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.embed import components
import vdj.io
import vdj.stats
bokeh.plotting.output_file('./point_endogenous_comparison.html')
# bokeh.io.output_notebook()


# Start by trying to figure out the details of picking the point mutants that
# make up the reference
endog_seqs = vdj.io.endogenous_seqs()
reference = endog_seqs['reference']
nt_idx = vdj.io.nucleotide_idx()

# Generate an empty dataframe and a color list from 5' to 3'
mut_df = pd.DataFrame([])
colors = bokeh.palettes.Category10_10
for key, val in endog_seqs.items():
    color_idx = 0   
    for i, b in enumerate(reference[1]):
        if b != val[1][i]:
            # Find what bases are present in the endogenous that aren't in ref
            base_letter = val[0][i]
            base_idx = nt_idx[val[0][i]] 

            # Define the name. 
            if i < 7:
                loc = 'Hept'
                pos = i + 1

            elif (i >= 7)  & (i < 19):
                loc = 'Spac'
                pos = (i + 1) - 7

            else:
                loc = 'Non'
                pos = (i + 1) - 19
            name = f'12{loc}{reference[0][i]}{pos}{base_letter}'
            color = colors[color_idx]

            # Assemble the mutant dataframe
            mut_info = {'position':i, 'base':base_letter, 'base_idx':base_idx,
                        'point_mutant': name, 'mutant':key, 'color':color, 'display_color':color}
            color_idx += 1
        else:
            mut_info = {'position':i, 'base':reference[0][i], 'base_idx':nt_idx[b],
                        'point_mutant':'No Mutation', 'mutant':key, 
                        'color':'#c2c2c2'}
        if (key != 'WT12rss') & (key != 'reference'):
            mut_df = mut_df.append(mut_info, ignore_index=True)
seq_source = ColumnDataSource(mut_df)

# %%
# Load the data sets
loops = pd.read_csv('../../../data/compiled_looping_events.csv')
dwell_all_data = pd.read_csv('../../../data/compiled_dwell_times.csv')
post_data = pd.read_csv('../../../data/pooled_cutting_probability_posteriors.csv')
post_data['color'] = 'slategrey'

# Restrict the posterior distributions.
post_data = post_data[(post_data['hmgb1']==80) & (post_data['salt']=='Mg')]

dfs = []
for g, d in post_data.groupby('mutant'):
    # Determine the class of mutant.
    if ('Spac' in g) | ('Hept' in g) | ('Non' in g):
        parsed_seq = vdj.io.mutation_parser(g) 
        mut_class = 'point'
    else:
        mut_class = 'endogenous'
    d = d.copy()
    d['class'] = mut_class
    d['alpha'] = 1
    dfs.append(d)
post_data = pd.concat(dfs)

# Split into endogenous and point mutations. 
post_endog = post_data[post_data['class']=='endogenous'] 
post_point = post_data[post_data['class']=='point'] 

# Rename the x and y such that I can loop effectively. 
post_point.rename(columns={'probability':'x_val', 'posterior':'y_val'}, 
                    inplace=True)

# Process the looping data into statistics
loops = loops[(loops['hmgb1']==80) & (loops['salt']=='Mg')]
pooled_df = pd.DataFrame()
rep_df = pd.DataFrame()
for g, d in loops.groupby('mutant'):
    # Determine the class of mutant. 
    if ('Spac' in g) | ('Hept' in g) | ('Non' in g):
        parsed_seq = vdj.io.mutation_parser(g) 
        mut_class = 'point'
        position = np.where(vdj.io.mutation_parser(g)['seq_idx'] != vdj.io.endogenous_seqs()['reference'][1])[0]

    else:
        mut_class = 'endogenous'
        position = [29]
    if len(position) == 1:
        pooled_df = pooled_df.append({'mutant': g, 'class':mut_class,
            'y_val':d['n_loops'].sum() / len(d['bead_idx']),
            'n_beads':len(d['bead_idx']), 'n_loops':d['n_loops'].sum(),
            'x_val':position[0], 'color': 'slategrey', 'alpha':1},
            ignore_index=True)
        for _g, _d in d.groupby(['replicate']):
            rep_df = rep_df.append({'mutant':g, 'class':mut_class,
            'y_val':_d['n_loops'].sum() / len(_d['bead_idx']),
            'n_beads':len(_d['bead_idx']), 'n_loops':_d['n_loops'].sum(),
            'x_val': position[0], 'color': 'slategrey', 'alpha':1},
            ignore_index=True)


#
#%% Process all datasets
# Keep only the HMGB1 = 80 mM and Mg
dwell_all_data = dwell_all_data[(dwell_all_data['hmgb1']==80) & 
                                (dwell_all_data['salt']=='Mg')]

# Separate into cut, unloop
dwell_cut_data = dwell_all_data[dwell_all_data['cut']==1]
dwell_unloop_data = dwell_all_data[dwell_all_data['cut']==0]

# Iterate through the dwell time data  and compute the ECDFS. 
dfs = []
for source in [dwell_all_data, dwell_cut_data, dwell_unloop_data]: 
    bin_dfs = []
    for g, d in source.groupby('mutant'):
        _df = pd.DataFrame()
        # Determine the mutant
        x, y = np.sort(d['dwell_time_min'].values), np.arange(0, len(d), 1) / len(d)
        y[-1] = 1 
        # Stagger the results so I can recreate a step plot
        staircase_y = np.empty(2 * len(d)) 
        staircase_x = np.empty(2 * len(d)) 
        staircase_y[0] = 0
        staircase_y[1::2] = y
        staircase_y[2::2] = y[:-1]
        staircase_x[::2] = x
        staircase_x[1::2] = x

        # Generate another point array so points can be plotted on arrays. 
        point_x = np.zeros(2 * len(d))
        point_y = np.zeros(2 * len(d))
        point_x[1::2] = x
        point_y[1::2] = y
        point_x[point_x == 0] = -1
        point_y[point_y == 0] = -1

        # Assemble the dataframe. 
        _df['x_val'] = staircase_x
        _df['point_x'] = point_x
        _df['y_val'] = staircase_y
        _df['point_y'] = point_y
        _df['mutant'] = g
        _df['color'] = 'slategrey'
        _df['alpha'] = 1;
        _df['level'] = 'underlay'
        if ('Spac' in g) | ('Hept' in g) | ('Non' in g):
            _df['class'] = 'point'
        else:
            _df['class'] = 'endogenous'
        bin_dfs.append(_df)     
    dwell_hist = pd.concat(bin_dfs)
    dfs.append(dwell_hist)
dwell_dist, cut_dist, unloop_dist = dfs

# Separate into point  and endogenous dists.  
dwell_all_point = ColumnDataSource(dwell_dist[dwell_dist['class']=='point']) 
dwell_all_endog  = ColumnDataSource(dwell_dist[dwell_dist['class']=='endogenous']) 
dwell_cut_point = ColumnDataSource(cut_dist[cut_dist['class']=='point'])
dwell_cut_endog = ColumnDataSource(cut_dist[cut_dist['class']=='endogenous'])
dwell_unloop_point = ColumnDataSource(unloop_dist[unloop_dist['class']== 'point'])
dwell_unloop_endog = ColumnDataSource(unloop_dist[unloop_dist['class']== 'endogenous'])

# Assemble teh posterior distribution cds
post_endog = ColumnDataSource(post_endog)
post_point = ColumnDataSource(post_point)
post_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})
leg_source = ColumnDataSource({'xs': [], 'ys': [], 'c':[], 'mutant':[], 'alpha':[]})

# Make blank dwell time cds for plotting. 
dwell_all_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})
dwell_cut_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})
dwell_unloop_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})

pooled_point = ColumnDataSource(pooled_df[pooled_df['class']=='point'])
rep_point = ColumnDataSource(rep_df[rep_df['class']=='point'])
pooled_endog = ColumnDataSource(pooled_df[pooled_df['class']=='endogenous'])
rep_endog = ColumnDataSource(rep_df[rep_df['class']=='endogenous'])


#%%

# Define filters
endog_filter = GroupFilter(column_name='mutant', group='')
target_filter = GroupFilter(column_name='mutant', group='')
dwell_all_filter = IndexFilter(indices=[])
dwell_cut_filter = IndexFilter(indices=[])
dwell_unloop_filter = IndexFilter(indices=[])
rep_filter = IndexFilter(indices=[])
pooled_filter = IndexFilter(indices=[])

# Define the many, many, (many) views
seq_view = CDSView(source=seq_source, filters=[endog_filter])
dwell_all_endog_view = CDSView(source=dwell_all_endog, filters=[endog_filter])
dwell_all_point_view = CDSView(source=dwell_all_point, filters=[dwell_all_filter])

dwell_cut_endog_view = CDSView(source=dwell_cut_endog, filters=[endog_filter])
dwell_cut_point_view = CDSView(source=dwell_cut_point, filters=[dwell_cut_filter])


dwell_unloop_endog_view = CDSView(source=dwell_unloop_endog, filters=[endog_filter])
dwell_unloop_point_view = CDSView(source=dwell_unloop_point, filters=[dwell_unloop_filter])


rep_point_view = CDSView(source=rep_point, filters=[rep_filter])
rep_endog_view = CDSView(source=rep_endog, filters=[endog_filter])

pooled_point_view = CDSView(source=pooled_point, filters=[pooled_filter])
pooled_endog_view = CDSView(source=pooled_endog, filters=[endog_filter])


post_endog_view = CDSView(source=post_endog, filters=[endog_filter])
post_point_view = CDSView(source=post_point)

# Define the dropdown for the interactivity
menu_dict = {m:m for m in mut_df['mutant'].unique()}
menu_dict['DFL161'] = "DFL16.1-5"
menu_dict['DFL1613'] = "DFL16.1-3"
menu = [(v,k) for k, v in menu_dict.items()]
sel = Dropdown(value='', label='Select Endogenous Sequence', menu=menu)

seq_hover_cb = """
var hover_mut = seq_source.data['point_mutant'][cb_data.index['1d'].indices[0]];
"""

sel_cb = """
var hover_mut = 'None';
"""

# Load the callback cod
with open('point_endogenous_comparison.js', 'r') as file:
    unified_code = file.read()
    sel_js = sel_cb + unified_code
    hover_js  = seq_hover_cb + unified_code

js_cbs = []
for cb in [sel_js, hover_js]:
    js_cbs.append(CustomJS(args={'endog_filter':endog_filter, 'leg_source':leg_source,
                    'dwell_all_filter':dwell_all_filter, 'dwell_cut_filter':dwell_cut_filter,
                    'dwell_unloop_filter':dwell_unloop_filter,
                    'rep_filter':rep_filter, 'pooled_filter':pooled_filter,
                    'seq_source':seq_source,  'seq_view':seq_view,
                    'dwell_all_endog':dwell_all_endog, 'dwell_all_endog_view': dwell_all_endog_view,
                    'dwell_cut_endog':dwell_cut_endog, 'dwell_cut_endog_view': dwell_cut_endog_view,
                    'dwell_unloop_endog':dwell_unloop_endog, 'dwell_unloop_endog_view': dwell_unloop_endog_view,
                    'dwell_all_point':dwell_all_point, 'dwell_all_point_view': dwell_all_point_view,
                    'dwell_cut_point':dwell_cut_point, 'dwell_cut_point_view': dwell_cut_point_view,
                    'dwell_unloop_point':dwell_unloop_point, 'dwell_unloop_point_view': dwell_unloop_point_view,
                    'dwell_all_blank': dwell_all_blank, 'dwell_cut_blank':dwell_cut_blank,
                    'dwell_unloop_blank':dwell_unloop_blank,
                    'rep_point':rep_point, 'rep_point_view':rep_point_view, 
                    'rep_endog':rep_endog, 'rep_endog_view':rep_endog_view,
                    'pooled_point':pooled_point, 'pooled_point_view':pooled_point_view, 
                    'pooled_endog':pooled_endog, 'pooled_endog_view':pooled_endog_view,
                    'post_endog': post_endog, 'post_endog_view':post_endog_view,
                    'post_point':post_point, 'post_view':post_point_view, 
                    'post_endog': post_endog, 'post_endog_view':post_endog_view,
                    'post_blank':post_blank,
                    'endog_sel':sel},
        code=cb))


# Define the figure canvases
seq_ax = bokeh.plotting.figure(width=600, height=50, x_range=[0, 30],
                                y_range=[-0.01, 0.1], tools=[''],
                                toolbar_location=None)
dwell_all_ax = bokeh.plotting.figure(width=400, height=250, x_axis_type='log', 
                                    x_range=[0.5, 80], y_range=[-0.05, 1.05],
                                    x_axis_label='paired complex dwell time [min]',
                                    y_axis_label = 'ECDF', title='all PCs',
                                    tools=[''])
dwell_unloop_ax = bokeh.plotting.figure(width=400, height=250, x_axis_type='log', 
                                    x_range=[0.5, 80], y_range=[0, 1],
                                    x_axis_label='paired complex dwell time [min]',
                                    y_axis_label = 'ECDF', title='unlooped PCs',
                                    tools=[''])
                            
dwell_cut_ax = bokeh.plotting.figure(width=400, height=250, x_axis_type='log', 
                                    x_range=[0.5, 80], y_range=[0, 1],
                                    x_axis_label='paired complex dwell time [min]',
                                    y_axis_label = 'ECDF', title='cleaved PCs',
                                    tools=[''])
loop_freq_ax = bokeh.plotting.figure(width=550, height=350,
                                    x_axis_label='reference nucleotide', 
                                    y_axis_label = 'paired complexes per bead',
                                    x_range=[0, 32], y_range=[0, 1],
                                    tools=[''], toolbar_location=None)
pcut_ax = bokeh.plotting.figure(width=550, height=350, x_axis_label='probability',
                                y_axis_label='posterior probability',
                                    title = 'PC cleavage probability',
                                    x_range=[0, 1],
                                tools=[''], toolbar_location=None)

# Define a legend axis and blank it 
leg_ax = bokeh.plotting.figure(width = 150, height= 700, tools=[''],
toolbar_location=None, x_range=[0, 1], y_range=[0, 1]) 

leg_ax.multi_line('xs', 'ys', color='c', line_width=10, legend='mutant', 
                  source=leg_source, alpha='alpha')
leg_ax.legend.location = 'center'
leg_ax.legend.background_fill_color = None
leg_ax.legend.label_text_font_size = '10pt'

# Format the sequence axis to not be colorful.
for ax in [seq_ax, leg_ax]:
    ax.xaxis.visible = False
    ax.yaxis.visible = False
    ax.background_fill_color = None
    ax.outline_line_color = None

# Set the ticker for the x axis o the sequence
ticks = np.arange(1, 30, 1)
ticks[-1] += 2
loop_freq_ax.ray(31, 0, angle=np.pi/2,length=20, line_color='white', line_width=15, alpha=0.5)
loop_freq_ax.xaxis.ticker = ticks
renamed_ticks = {int(t):s for t, s in zip(ticks, list(reference[0]))}
renamed_ticks[31] = 'endo'
loop_freq_ax.xaxis.major_label_overrides = renamed_ticks

# Add the variant sequences
variant_seq = bokeh.models.glyphs.Text(x='position', y=0, text='base',
text_color='color', text_font='Courier', text_font_size='26pt')
sequence = seq_ax.add_glyph(seq_source, variant_seq, view=seq_view)

# Add a hover interaction to the seq ax
loop_freq_ax.triangle(31, 'y_val', source=rep_endog, 
                      view=rep_endog_view, color='slategrey', alpha=0.75, 
                      size=8)
loop_freq_ax.circle(31, 'y_val', source=pooled_endog, 
                      view=pooled_endog_view,
                      fill_color='slategrey', line_color='black',  
                      size=10, line_width=1)

# Plot the ECDFS of the endogenous samples
dwell_all_rend = dwell_all_ax.line('x_val', 'y_val', line_width=2, color='slategrey', 
                   source=dwell_all_endog, view=dwell_all_endog_view,
                   level='overlay')
dwell_cut_rend = dwell_cut_ax.line('x_val', 'y_val', line_width=2, color='slategrey', 
                   source=dwell_cut_endog, view=dwell_cut_endog_view)
dwell_unloop_rend = dwell_unloop_ax.line('x_val', 'y_val', line_width=2, color='slategrey', 
                   source=dwell_unloop_endog, view=dwell_unloop_endog_view)

loop_freq_ax.triangle('x_val', 'y_val', source=rep_point, 
                      view=rep_point_view, color='color', alpha='alpha', 
                      size=8)
loop_rend = loop_freq_ax.circle('x_val', 'y_val', source=pooled_point, 
                      view=pooled_point_view,
                      color='color', line_color='black',  
                      size=10, line_width=1,
                      alpha='alpha')

# Plot the endogenous posterior distributions
pcut_ax.line('probability', 'posterior', source=post_endog, view=post_endog_view,
            line_width=2, color='slategrey', alpha=0.75)

pcut_ax.varea('probability', y1=0, y2='posterior', source=post_endog, view=post_endog_view,
            color='slategrey', alpha=0.25)

# # Plot the dwell distribution for the point mutants. 
dwell_all_ax.multi_line('xs', 'ys', source=dwell_all_blank, color='c',
                        line_width=2, alpha='alpha') 
dwell_cut_ax.multi_line('xs', 'ys', source=dwell_cut_blank, color='c',
                        line_width=2, alpha='alpha') 
dwell_unloop_ax.multi_line('xs', 'ys', source=dwell_unloop_blank, color='c',
                            line_width=2, alpha='alpha')
post_rend = pcut_ax.multi_line('xs', 'ys', source=post_blank, line_width=2, color='c',
                    alpha='alpha')


hover = HoverTool(renderers=[sequence], callback=js_cbs[1],
tooltips=[('mutation', '@point_mutant')])

seq_ax.add_tools(hover)
sel.js_on_change('value', js_cbs[0])


spacer = Div(text='<br/>')
sel_row = bokeh.layouts.row(sel, seq_ax, sizing_mode='scale_width')
loop_col = bokeh.layouts.column(spacer, loop_freq_ax, spacer, pcut_ax)
dwell_col = bokeh.layouts.column(dwell_unloop_ax, dwell_cut_ax, dwell_all_ax)
complete_row = bokeh.layouts.row(leg_ax, loop_col, dwell_col)

lay = bokeh.layouts.column(sel_row, complete_row)

# Set the theme. 
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


# FOrmat legend details. 
loop_freq_ax.legend.click_policy = 'hide'
# loop_freq_ax.legend.title_text = 'click to hide'
bokeh.io.save(lay)



#%%
