#%%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.layouts
import bokeh.plotting
import matplotlib.patches as patches
import vdj.io
import vdj.viz
vdj.viz.plotting_style_bokeh()
bokeh.io.output_notebook()

# Upload data on coding flank-relevant RSSs
cf_muts = ['WT12rss', '12CodC6A', '12SpacC1A', 'V4-55']

loop = pd.read_csv('../../data/compiled_loop_freq_bs.csv')
loop = loop[(loop['mutant'].isin(cf_muts)) & (loop['hmgb1']==80) & (loop['salt']=='Mg')]
loop = loop.replace(to_replace='WT12rss', value='V4-57-1 (ref)')

dwell = pd.read_csv('../../data/compiled_dwell_times.csv')
dwell = dwell[(dwell['mutant'].isin(cf_muts)) & (dwell['hmgb1']==80) & (dwell['salt']=='Mg')]
dwell = dwell.replace(to_replace='WT12rss', value='V4-57-1 (ref)')
dwell_cut = dwell[dwell['cut']==1].copy()
dwell_unloop = dwell[dwell['cut']==0].copy()

cut_posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv')
cut_posts = cut_posts[(cut_posts['mutant'].isin(cf_muts)) & (cut_posts['hmgb1']==80) & (cut_posts['salt']=='Mg')]
cut_posts = cut_posts.replace(to_replace='WT12rss', value='V4-57-1 (ref)')

#%%
# Create dwell ECDFs
dfs = []
for source in [dwell, dwell_cut, dwell_unloop]:
    dist_df = []
    for g,d in source.groupby('mutant'):
        x,y = np.sort(d['dwell_time_min'].values), np.arange(0, len(d), 1) / len(d)
        y[-1] = 1
        _df = pd.DataFrame()
        _df['x'] = x
        _df['y'] = y
        _df['mutant'] = g
        dist_df.append(_df)
    dwell_dist = pd.concat(dist_df)
    dfs.append(dwell_dist)
dwell_dist, cut_dist, unloop_dist = dfs

#%%
ax_conf_int = bokeh.plotting.figure(height=100, width=570)
ax_loop = bokeh.plotting.figure(height=200, width=570, 
                            x_range=[0,0.5], y_range=[0.1,0.65],
                            x_axis_label='looping frequency')
ax_loop.yaxis.visible=False

ax_dwell_cut = bokeh.plotting.figure(height=190, width=210, x_axis_type='log',
                                    x_range=[0.5,80], y_axis_label='ECDFs',
                                    x_axis_label='dwell time [min]')
ax_dwell_unloop = bokeh.plotting.figure(height=190, width=180, x_axis_type='log',
                                x_range=[0.5, 80], x_axis_label='dwell time [min]')
ax_dwell_all = bokeh.plotting.figure(height=190, width=180, x_axis_type='log',
                                x_range=[0.5, 80], x_axis_label='dwell time [min]')

ax_posts = bokeh.plotting.figure(height=200, width=570, 
                                x_range=[0, 1.0], x_axis_label='cutting probability',
                                y_axis_label='posterior')

seqs = {'V4-57-1 (ref)':0.25, '12CodC6A':0.5}
colors = {'V4-57-1 (ref)':'slategrey', '12CodC6A':'dodgerblue'}

for seq in seqs:
    ax_loop.triangle(loop[loop['mutant']==seq]['loops_per_bead'], seqs[seq], 
                    size=15, fill_color='white', line_color=colors[seq], level='overlay')
    for g,d in loop[loop['mutant']==seq].groupby('percentile'):
        ax_loop.rect((d['high'] + d['low']) / 2, seqs[seq], d['high'] - d['low'], 
                    0.15, color=colors[seq], alpha=0.15)
    ax_dwell_cut.line(cut_dist[cut_dist['mutant']==seq]['x'], 
                        cut_dist[cut_dist['mutant']==seq]['y'],
                        line_width=2, line_color=colors[seq])
    ax_dwell_unloop.line(unloop_dist[unloop_dist['mutant']==seq]['x'], 
                        unloop_dist[unloop_dist['mutant']==seq]['y'],
                        line_width=2, line_color=colors[seq])
    ax_dwell_all.line(dwell_dist[dwell_dist['mutant']==seq]['x'], 
                        dwell_dist[dwell_dist['mutant']==seq]['y'],
                        line_width=2, line_color=colors[seq])
    ax_posts.line(cut_posts[cut_posts['mutant']==seq]['probability'],
                cut_posts[cut_posts['mutant']==seq]['posterior'],
                line_width=2, line_color=colors[seq])
    ax_posts.varea(cut_posts[cut_posts['mutant']==seq]['probability'],
                0, cut_posts[cut_posts['mutant']==seq]['posterior'],
                color=colors[seq], alpha=0.4)

dwell_row = bokeh.layouts.row(ax_dwell_cut, ax_dwell_unloop, ax_dwell_all)
column = bokeh.layouts.column(ax_loop, dwell_row, ax_posts)
bokeh.io.show(column)

#%%
