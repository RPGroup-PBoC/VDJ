import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, 
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool,
                          LinearColorMapper)
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.embed import components
import vdj.io
import vdj.stats
from bokeh.transform import transform
import bokeh.palettes
bokeh.plotting.output_file('./point_voyager.html')

# Load the necessary data sets
dwell_times = pd.read_csv('../../../data/compiled_dwell_times.csv')
posteriors = pd.read_csv('../../../data/pooled_cutting_probability_posteriors.csv')
loops = pd.read_csv('../../../data/compiled_looping_events.csv')

# Keep the data from the reference sequence
dwell_ref = dwell_times[(dwell_times['mutant']=='WT12rss') & 
                       (dwell_times['hmgb1']==80) & (dwell_times['salt']=='Mg')]
cut_ref = dwell_ref[dwell_ref['cut']==1]
post_ref = posteriors[(posteriors['mutant']=='WT12rss') & 
                             (posteriors['hmgb1']==80) & (posteriors['salt']=='Mg')]
loops_ref = loops[(loops['mutant']=='WT12rss') & 
                    (loops['hmgb1']==80) & (loops['salt']=='Mg')]

# Identify the endogenous muts
dfs  = []
for i, df in enumerate([dwell_times, posteriors, loops]):
    for g, d in df.groupby(['mutant']):
        if ('12Spac' not in g) & ('12Hept' not in g) & ('12Non' not in g):
            mut_class = 'endogenous'
        else:
            mut_class = 'point' 
       
        df.loc[df['mutant']==g, 'class'] = mut_class
    df = df[(df['class']=='point') & (df['salt']=='Mg') & (df['hmgb1']==80)]
    dfs.append(df)
dwell_times, posteriors, loops = dfs
cut_dwells = dwell_times[dwell_times['cut']==1]

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

# Do the seame for the reference
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
rep_ref['y] = np.random.normal(0,0.05, len(rep_df))

# ##############################################################################
# GENERATE HISTOGRAMS OF DWELL TIMES 
# ##############################################################################
# Generate the histogrammed dwell times.
# bins = np.linspace(0, dwell_times['dwell_time_min'].max(), 25)
# dfs = []
# for source in [dwell_times, cut_dwells]: 
#     bin_dfs = []
#     for g, d in source.groupby('mutant'):
#         print('AAHHHHHHH')
#         hist, bins = np.histogram(d['dwell_time_min'], bins=bins)
#         _df = pd.DataFrame()
#         _df['top'] = hist
#         _df['bottom'] = 0
#         _df['left'] = bins[1:]
#         _df['right'] = bins[:-1]
#         _df['mutant'] = g
#         bin_dfs.append(_df)
#     dwell_hist = pd.concat(bin_dfs)
#     dfs.append(dwell_hist)
# dwell_hist, cut_hist = dfs

# # Do the same for the reference sequence
# dfs = []
# for source in [dwell_ref, cut_ref]: 
#     bin_dfs = []
#     for g, d in source.groupby('mutant'):
#         hist, bins = np.histogram(d['dwell_time_min'], bins=bins)
#         _df = pd.DataFrame()
#         _df['top'] = hist
#         _df['bottom'] = 0
#         _df['left'] = bins[1:]
#         _df['right'] = bins[:-1]
#         _df['mutant'] = g
#         bin_dfs.append(_df)
#     dwell_hist = pd.concat(bin_dfs)
#     dfs.append(dwell_hist)
# dwell_hist_ref, cut_hist_ref = dfs

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
                'pos': loc + 1,
                'base_idx': nt_idx[base],
                'base': base,
                'diff': diff,
                'med_dwell':med_dwell}, ignore_index=True)
dwell_source = ColumnDataSource(dwell_mat)

# ##############################################################################
# COMPUTE THE DIFFERENCE IN LOOP FREQUENCY
# ##############################################################################
loops_ref = pooled_ref['loops_per_bead']
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
                'pos': loc + 1,
                'base_idx': nt_idx[base],
                'base': base,
                'diff': diff,
                'med_dwell':med_dwell}, ignore_index=True)
d


# Set up a figure. 
p = bokeh.plotting.figure(height=100)
palette = bokeh.palettes.RdBu11
dwell_color = LinearColorMapper(palette=palette, low=-3, high=3)
p.rect(np.arange(0, 28, 1), ref_seq, width=1, height=1, fill_color='white',
        line_color='black')
p.circle(np.arange(0, 28, 1), ref_seq, line_color='slategrey')

p.rect('pos', 'base_idx', width=1, height=1, source=dwell_source, 
        fill_color=transform('diff', dwell_color))

bokeh.io.show(p)