#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import vdj.io 
import vdj.viz
vdj.viz.plotting_style()

# Load the data with long-form looping events and restrict to relevant sets.
data = pd.read_csv('../../data/compiled_looping_events.csv')
data = data[(data['salt']=='Mg') & (data['hmgb1']==80)]

# Load the dwell time data
dwell = pd.read_csv('../../data/compiled_dwell_times.csv')
dwell = dwell[(dwell['salt']=='Mg') & (dwell['hmgb1']==80)]

# Load all cutting probability estimates taking gaussian approximation.
cut_data = pd.read_csv('../../data/pooled_cutting_probability.csv')
cut_data = cut_data[(cut_data['hmgb1'] == 80) & (cut_data['salt']=='Mg')]

# Compute the number of loops per bead
counts = data.groupby(['mutant'])[['n_loops']].agg(('sum', 'count')).reset_index()
counts['loops_per_bead'] = counts['n_loops']['sum'] / counts['n_loops']['count']

#Compute the median dwell time
median_dwell = dwell.groupby('mutant')['dwell_time_min'].median().reset_index()

# Load the reference sequence
ref = vdj.io.endogenous_seqs()
ref_seq = ref['WT12rss'][0]
ref_idx = ref['WT12rss'][1]

# Compute the wild-type cutting probability and loop freq
wt_pcut = cut_data[cut_data['mutant']=='WT12rss']['mode'].values[0]
wt_freq = counts[counts['mutant']=='WT12rss']['loops_per_bead'].values[0]
wt_dwell = median_dwell[median_dwell['mutant']=='WT12rss']['dwell_time_min'].values[0]

# Compute the relative difference. 
cut_data['rel_diff'] = cut_data['mode'] - wt_pcut
counts['rel_diff'] = counts['loops_per_bead'] - wt_freq
median_dwell['rel_diff'] = median_dwell['dwell_time_min'] - wt_dwell

#%%
# Make a new data frame and connect the two statistics. 
df = pd.DataFrame([])
for g, d in cut_data.groupby('mutant'):
    # Get the looping frequencies for the same mutant. 
    freq = counts[counts['mutant']==g].copy()
    _dwell = median_dwell[median_dwell['mutant']==g].copy()

    # Parse the info about the mutant. 
    parsed = vdj.io.mutation_parser(g)
    if parsed['n_muts'] == 1:
        exp = 'point'

        # Get the mutated base and position.
        pos = np.argmax((parsed['seq_idx'] != ref_idx).astype(int))
        mut = parsed['seq'][pos]

        # Determine the location of the mutation.
        if 'Hept' in g:
            region = 'Heptamer'
        elif 'Spac' in g:
            region = 'Spacer'
        else:
            region = 'Nonamer'
    elif parsed['n_muts'] == 0:
        exp = 'reference'
        pos = -1
        mut = 'X'
        region = 'X'
    else:
        exp = 'endogenous'
        pos = -1
        mut = 'X'
        region = 'X'

    # Generate the new data frame. 
    df = df.append({'freq':freq['loops_per_bead'].values[0],
                   'p_cut':d['mode'].values[0],
                   'dwell':_dwell['dwell_time_min'].values[0],
                   'freq_diff':freq['rel_diff'].values[0],
                   'prob_diff':d['rel_diff'].values[0],
                   'dwell_diff':_dwell['rel_diff'].values[0],
                   'mutant': g,
                   'kind': exp,
                   'region':region,
                   'position':pos + 1,
                   'mutation': mut}, ignore_index=True)
                   
#%%
# Set up the figure canvas
fig, ax = plt.subplots(3, 3, figsize=(5, 5))

# Turn off unnecessary axes
ax[0, 0].axis('off')
ax[0, 1].axis('off')
ax[0, 2].axis('off')
ax[1, 1].axis('off')
ax[1, 2].axis('off')
ax[2, 2].axis('off')

# Turn off the ticks
ax[1, 0].set_xticks([])
ax[2, 1].set_yticks([])

# Plot the dashed zero lines
populated_ax = [ax[1, 0]

#%%
# loop_marg = fig.add_subplot(gs[1:5, -1])
_marg.set_xticks([])
cut_marg.set_yticks([])
loop_marg.set_xticks([])
loop_marg.set_yticks([])

# Set the zero lines. 
ax.hlines(0, -0.5, 0.5, 'k', linestyle='--')
ax.vlines(0, -0.7, 0.7, 'k', linestyle='--')

# Change the limits
ax.set_xlim([-0.5, 0.5])
ax.set_ylim([-0.7, 0.7])
ax.set_ylabel('change in\ncutting probability')
ax.set_xlabel('change in\nlooping frequency')

# Plot the point mutations as open circles
points = df[df['kind']=='point']

# Assign colors for regions
region_colors = ['rebeccapurple', 'dodgerblue', 'tomato']

i = 0
zorder = 4 
for g, d in points.groupby('region'):
    # Compute the histograms of each for the marginal distributions. 
    freq_bins= np.linspace(-0.5, 0.5, 15)
    prob_bins = np.linspace(-0.7, 0.7, 15)
    # freq_marg.hist(d['freq_diff'], bins=20, density=True, color=region_colors[i],
                #    alpha=0.5)
    cut_marg.hist(d['prob_diff'], bins=20, color=region_colors[i],
                   alpha=0.5, orientation='horizontal')
    ax.plot(d['freq_diff'], d['prob_diff'], 'o', ms=4, 
            color=region_colors[i], alpha=0.78, label=g, zorder=zorder)
    i += 1
    zorder -= 1 
    # ax.annotate(g['mutation'], xy=(d['freq_diff'], d['prob_diff']))
# ax.legend()
plt.savefig('./cutting_vs_looping_diff.pdf', bbox_inches='tight')
#%%


#%%
