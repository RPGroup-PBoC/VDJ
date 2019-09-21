#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

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
fig = plt.figure(constrained_layout=False, figsize=(9,6))
fig.xticklabels([])
fig.set_yticklabels([])
gs_ci = fig.add_gridspec(nrows=1, ncols=3)
ax_conf_int = fig.add_subplot(gs_ci[:, :])

gs = fig.add_gridspec(2, 3)
ax_loop = fig.add_subplot(gs[0,0])
ax_loop.set_xticklabels([])
ax_loop.set_xlim([0.2, 0.7])
ax_loop.set_ylim([0, 0.3])
ax_loop.set_yticklabels([0, 0.1, 0.2, 0.3])

ax_posts = fig.add_subplot(gs[0,1:])

ax_dwell_cut = fig.add_subplot(gs[1,0])
ax_dwell_cut.set_xscale('log')

ax_dwell_unloop = fig.add_subplot(gs[1,1])
ax_dwell_unloop.set_xscale('log')

ax_dwell_all = fig.add_subplot(gs[1,2])
ax_dwell_all.set_xscale('log')

seqs = {'V4-57-1 (ref)':0.25, '12CodC6A':0.5}
colors = {'V4-57-1 (ref)':'slategrey', '12CodC6A':'dodgerblue'}

for seq in seqs:
    ax_loop.scatter(seqs[seq], loop[loop['mutant']==seq]['loops_per_bead'].values[0], 
                    s=15, color='white', edgecolor=colors[seq], zorder=10)
    for g,d in loop[loop['mutant']==seq].groupby('percentile'):
        rect = patches.Rectangle([seqs[seq], d['low']], 0.15,
                                d['high'].values[0] - d['low'].values[0],
                                color=colors[seq], alpha=0.2)
        ax_loop.add_patch(rect)
    ax_dwell_cut.plot(cut_dist[cut_dist['mutant']==seq]['x'], 
                        cut_dist[cut_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_dwell_unloop.plot(unloop_dist[unloop_dist['mutant']==seq]['x'], 
                        unloop_dist[unloop_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_dwell_all.plot(dwell_dist[dwell_dist['mutant']==seq]['x'], 
                        dwell_dist[dwell_dist['mutant']==seq]['y'],
                        lw=2, color=colors[seq])
    ax_posts.plot(cut_posts[cut_posts['mutant']==seq]['probability'],
                    cut_posts[cut_posts['mutant']==seq]['posterior'],
                    lw=2, color=colors[seq])
    ax_posts.fill_between(cut_posts[cut_posts['mutant']==seq]['probability'],
                            0, cut_posts[cut_posts['mutant']==seq]['posterior'],
                            color=colors[seq], alpha=0.4)

#%%
ax_conf_int = bokeh.plotting.figure(height=50, width=570, 
                                    x_range=[0.125, 0.805], y_range=[0, 0.1])
ax_conf_int.xaxis.visible = False
ax_conf_int.yaxis.visible = False
ax_conf_int.background_fill_color = 'white'

ax_loop = bokeh.plotting.figure(height=200, width=210, 
                                x_range=[0.1,0.65], y_range=[0.0,0.35],
                                y_axis_label='looping frequency')
ax_loop.xaxis.visible=False
ax_loop.yaxis.ticker=[0, 0.1, 0.2, 0.3]

ax_dwell_cut = bokeh.plotting.figure(height=190, width=210, x_axis_type='log',
                                    x_range=[0.5,80], y_axis_label='ECDFs',
                                    x_axis_label='dwell time [min]')
ax_dwell_unloop = bokeh.plotting.figure(height=190, width=180, x_axis_type='log',
                                    x_range=[0.5, 80], x_axis_label='dwell time [min]')
ax_dwell_all = bokeh.plotting.figure(height=190, width=180, x_axis_type='log',
                                    x_range=[0.5, 80], x_axis_label='dwell time [min]')

ax_posts = bokeh.plotting.figure(height=200, width=360, 
                                x_range=[0, 1.0], x_axis_label='cutting probability',
                                y_axis_label='posterior')
ax_posts.yaxis.ticker = [0, 0.01, 0.02]

seqs = {'V4-55':0.25, '12SpacC1A':0.5}
colors = {'V4-55':'slategrey', '12SpacC1A':'dodgerblue'}

ci_widths = np.arange(0.03, 0.03 * (len(loop['percentile'].unique()) + 1), 0.03)
ci_centers = np.arange(0.015, 0.015 * (len(loop['percentile'].unique()) + 1), 0.015)
ci_text = np.arange(0.015, 0.015 + 0.030 * (len(loop['percentile'].unique())), 0.03)
ci_y = 0.08
ci_perc = np.sort(loop['percentile'].unique())
ci_df = pd.DataFrame({'percentile':ci_perc, 'y':ci_y, 'text':ci_text, 
                    'center':ci_centers, 'width':ci_widths})
ci_df['perc_text'] = list(str(int(perc)) for perc in ci_perc)
text_source = ColumnDataSource(ci_df)

for seq in seqs:
    ax_loop.triangle(seqs[seq], loop[loop['mutant']==seq]['loops_per_bead'], 
                    size=15, fill_color='white', line_color=colors[seq], level='overlay')
    for g,d in loop[loop['mutant']==seq].groupby('percentile'):
        ax_loop.rect(seqs[seq], (d['high'] + d['low']) / 2, 0.15,
                    d['high'] - d['low'], color=colors[seq], alpha=0.2)
        ax_conf_int.rect(seqs[seq] + ci_df[ci_df['percentile']==g]['center'], 0.05,
                        ci_df[ci_df['percentile']==g]['width'], 0.025, color=colors[seq],
                        alpha=0.2)
#        Text(x='ci_text', y='y', text='perc_text', source=text_source)
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

top_row = bokeh.layouts.row(ax_loop, ax_posts)
dwell_row = bokeh.layouts.row(ax_dwell_cut, ax_dwell_unloop, ax_dwell_all)
column = bokeh.layouts.column(ax_conf_int, top_row, dwell_row)
bokeh.io.show(column)

#%%


#%%
