#%% 
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, Text,
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool,
                          IndexFilter)
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

mut_df = pd.DataFrame([])
colors = bokeh.palettes.Paired11[1:]
for key, val in endog_seqs.items():
    # Find where the reference is different from the endog 
    color_idx = 0   
    for i, b in enumerate(reference[1]):
        if b != val[1][i]:
            # Find what bases are present in the endogenous that aren't in ref
            base_letter = val[0][i]
            base_idx = nt_idx[val[0][i]] 

            # Define the name. 
            if i < 8:
                loc = 'Hept'
                pos = i + 1
            elif (i >= 8)  & (i < 19):
                loc = 'Spac'
                pos = (i + 1) - 7
            else:
                loc = 'Non'
                pos = (i + 1) - 19
            name = f'12{loc}{reference[0][i]}{pos}{base_letter}'
            color = colors[color_idx]

            # Assemble the mutant dataframe
            mut_info = {'position':i, 'base':base_letter, 'base_idx':base_idx,
                        'mutant': name, 'endogenous':key, 'color':color}
            color_idx += 1
        else:
            mut_info = {'position':i, 'base':reference[0][i], 'base_idx':nt_idx[b],
                        'mutant':'No Mutation', 'endogenous':key, 
                        'color':'#c2c2c2'}
        mut_df = mut_df.append(mut_info, ignore_index=True)
seq_source = ColumnDataSource(mut_df)

# %%
# Load the data sets
loops = pd.read_csv('../../../data/compiled_looping_events.csv')
dwell_all_data = pd.read_csv('../../../data/compiled_dwell_times.csv')

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
            'loops_per_bead':d['n_loops'].sum() / len(d['bead_idx']),
            'n_beads':len(d['bead_idx']), 'n_loops':d['n_loops'].sum(),
            'position':position[0], 'color': 'slategrey'},
            ignore_index=True)
        for _g, _d in d.groupby(['replicate']):
            rep_df = rep_df.append({'mutant':g, 'class':mut_class,
            'loops_per_bead':_d['n_loops'].sum() / len(_d['bead_idx']),
            'n_beads':len(_d['bead_idx']), 'n_loops':_d['n_loops'].sum(),
            'position': position[0], 'color': 'slategrey'},
            ignore_index=True)
pooled_df['y'] = 0.5 
rep_df['y'] = np.random.normal(0.5,0.05, len(rep_df))
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
        _df['dwell_time'] = x
        _df['ecdf'] = y
        _df['mutant'] = g
        _df['color'] = 'slategrey'
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


# Make blank dwell time cds for plotting. 
dwell_all_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[]})
dwell_cut_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[]})
dwell_unloop_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[]})

pooled_point = ColumnDataSource(pooled_df[pooled_df['class']=='point'])
rep_point = ColumnDataSource(rep_df[rep_df['class']=='point'])
pooled_endog = ColumnDataSource(pooled_df[pooled_df['class']=='endogenous'])
rep_endog = ColumnDataSource(rep_df[rep_df['class']=='endogenous'])

#%%

# Define filters
endog_filter = GroupFilter(column_name='endogenous', group='')
target_filter = GroupFilter(column_name='mutant', group='')
dwell_all_filter = IndexFilter(indices=[])
dwell_cut_filter = IndexFilter(indices=[])
dwell_unloop_filter = IndexFilter(indices=[])
rep_filter = IndexFilter(indices=[])
pooled_filter = IndexFilter(indices=[])

# Define  the many, many views
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
pooled_endog_view = CDSView(source=pooled_endog, filters=[pooled_filter])

# Define the dropdown for the interactivity
sel = Select(options=list(mut_df['endogenous'].unique()))


# Load the callback code
with open('point_endogenous_comparison.js', 'r') as file:
    cb_code = file.read()

cb = CustomJS(args={'endog_filter':endog_filter, 'target_filter':target_filter, 
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
                    'endog_sel':sel},
    code=cb_code) 

sel.js_on_change('value', cb)


# Define the figure canvases
seq_ax = bokeh.plotting.figure(width=700, height=50, x_range=[0, 30],
                                tools=[''], toolbar_location=None,
                                y_range=[-0.01, 0.1])
                                
dwell_all_ax = bokeh.plotting.figure(width=600, height=200, x_axis_type='log', 
                                    x_range=[0.8, 80], 
                                    x_axis_label='paired complex dwell time [min]',
                                    y_axis_label = 'ECDF', title='all PCs')
dwell_unloop_ax = bokeh.plotting.figure(width=300, height=200, x_axis_type='log', 
                                    x_range=[0.8, 80], 
                                    x_axis_label='paired complex dwell time [min]',
                                    y_axis_label = 'ECDF', title='unlooped PCs')

dwell_cut_ax = bokeh.plotting.figure(width=300, height=200, x_axis_type='log', 
                                    x_range=[0.8, 80], 
                                    x_axis_label='paired complex dwell time [min]',
                                    y_axis_label = 'ECDF', title='cleaved PCs')
loop_freq_ax = bokeh.plotting.figure(width=700, height=400,
                                    x_axis_label='position of mutation', 
                                    y_axis_label = 'paired complexes per bead',
                                    x_range=[0, 30], y_range=[0, 0.8])
seq_ax.xaxis.visible = False
seq_ax.yaxis.visible = False
seq_ax.background_fill_color = None
seq_ax.outline_line_color = None

# Add the variant sequences
variant_seq = bokeh.models.glyphs.Text(x='position', y=0, text='base',
text_color='color', text_font='Courier', text_font_size='28pt')
seq_ax.add_glyph(seq_source, variant_seq, view=seq_view)

# Plot the ECDFS of the endogenous samples
dwell_all_ax.step('dwell_time', 'ecdf', line_width=1, color='slategrey', 
                   source=dwell_all_endog, view=dwell_all_endog_view)
dwell_all_ax.circle('dwell_time', 'ecdf', line_width=1, fill_color='white', color='slategrey', 
                   source=dwell_all_endog, view=dwell_all_endog_view)

dwell_cut_ax.step('dwell_time', 'ecdf', line_width=1, color='slategrey', 
                   source=dwell_cut_endog, view=dwell_cut_endog_view)
dwell_cut_ax.circle('dwell_time', 'ecdf', line_width=1, fill_color='white', color='slategrey', 
                   source=dwell_cut_endog, view=dwell_cut_endog_view)

dwell_unloop_ax.step('dwell_time', 'ecdf', line_width=1, color='slategrey', 
                   source=dwell_unloop_endog, view=dwell_unloop_endog_view)
dwell_unloop_ax.circle('dwell_time', 'ecdf', line_width=1, fill_color='white', color='slategrey', 
                   source=dwell_unloop_endog, view=dwell_unloop_endog_view)


# Plot the looping freqs for endogenous
loop_freq_ax.triangle('position', 'loops_per_bead', source=rep_endog, 
                      view=rep_endog_view, color='slategrey', alpha=0.75, 
                      size=8)
loop_freq_ax.circle('position', 'loops_per_bead', source=pooled_endog, 
                      view=pooled_endog_view,
                      line_color='slategrey',fill_color='white',  
                      size=10)

#
# # Plot the dwell distribution for the point mutants. 
dwell_all_ax.multi_line('xs', 'ys', source=dwell_all_blank, color='c', line_width=1)                   
dwell_all_ax.circle('xs', 'ys', line_width=1, 
                   source=dwell_all_blank,
                   color='c', alpha=0.5)

dwell_cut_ax.multi_line('xs', 'ys', source=dwell_cut_blank, color='c', line_width=1)                   
dwell_unloop_ax.multi_line('xs', 'ys', source=dwell_unloop_blank, color='c', line_width=1)                   


loop_freq_ax.triangle('position', 'loops_per_bead', source=rep_point, 
                      view=rep_point_view, color='color', alpha=0.5, 
                      size=8)
loop_freq_ax.circle('position', 'loops_per_bead', source=pooled_point, 
                      view=pooled_point_view,
                      line_color='color', fill_color='white',  
                      size=10)


# dwell_cut_ax.step('dwell_time', 'ecdf', line_width=1, 
#                    source=dwell_cut_point, view=dwell_cut_point_view,
#                    color='color')
# dwell_cut_ax.circle('dwell_time', 'ecdf', line_width=1, 
#                    source=dwell_cut_point, view=dwell_cut_point_view,
#                    color='color', alpha=0.5)
# dwell_unloop_ax.step('dwell_time', 'ecdf', line_width=1, 
#                    source=dwell_cut_point, view=dwell_cut_point_view,
#                    color='color')
# dwell_unloop_ax.circle('dwell_time', 'ecdf', line_width=1, 
#                    source=dwell_unloop_point, view=dwell_unloop_point_view,
#                    color='color', alpha=0.5)




dwell_row = bokeh.layouts.row(dwell_unloop_ax, dwell_cut_ax)
dwell_col = bokeh.layouts.column(dwell_row, dwell_all_ax)
seq_col = bokeh.layouts.column(sel, seq_ax, loop_freq_ax)
lay = bokeh.layouts.row(seq_col, dwell_col)

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
bokeh.io.save(lay)



#%%
