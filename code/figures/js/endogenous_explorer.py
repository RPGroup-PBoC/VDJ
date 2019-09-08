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
bokeh.io.output_notebook()  

#%%
# Load the necessary data sets
dwell_times = pd.read_csv('../../../data/compiled_dwell_times.csv')
# posts = pd.read_csv('../../../data/coutting_rate_posteriors')
# fates = pd.read_csv('../../../data/compiled_cutting_events.csv')


# Identify the endogenous muts
dfs  = []
for i, df in enumerate([dwell_times]):
    for g, d in df.groupby(['mutant']):
        if ('12Spac' not in g) & ('12Hept' not in g) & ('12Non' not in g):
            mut_class = 'endogenous'
        else:
            mut_class = 'point' 
       
        df.loc[df['mutant']==g, 'class'] = mut_class
    df = df[(df['class']=='endogenous') & (df['salt']=='Mg') & (df['hmgb1']==80)]
    dfs.append(df)
dwell_times = dfs[0]

#%%
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



# %% 
# ##############################################################################
# CALLBACK DEFINITION
# ##############################################################################
with open('endogenous_explorer.js') as f:
    cb = f.read()

#%%
# Set up the dropdown with the mutations
mut_sel = Dropdown(label = "Select Endogenous Mutation", 
            menu=[(m, m) for m in dwell_times['mutant'].unique()])

# Define a global mutant filter.
filter = GroupFilter(column_name="mutant", group="WT12rss")
view = CDSView(source=seq_variant, filters=[filter])

# Set up sequence axis
ax_seq = bokeh.plotting.figure(height=40, y_range=[0, 0.1])
ax_seq.xaxis.visible = False
ax_seq.yaxis.visible = False
ax_seq.grid.visible = False


seq_glyph = bokeh.models.glyphs.Text(x='pos', y='y', text='base', text_font='Courier',

                                    text_color='#95a3b2', text_font_size='28px')
ax_seq.add_glyph(invariant, seq_glyph)
var_glyph = bokeh.models.glyphs.Text(x='pos', y='y', text='base', text_font='Courier',
                                    text_color='tomato', text_font_size='28px')
ax_seq.add_glyph(seq_view, var_glyph)


ax_loop = bokeh.plotting.figure(height=100)
ax_dwell = bokeh.plotting.figure(height=250)
ax_cut = bokeh.plotting.figure(height=250)
# Set up the callback
callback = CustomJS(args={'drop':mut_sel, 'view':view, 'seqs':variant,
                          'filter':filter})
mut_sel.js_on_change('value', callback)


lay = bokeh.layouts.column(mut_sel, ax_seq) #, ax_loop, ax_dwell, ax_cut)
bokeh.io.show(lay)

#%%



#%%
