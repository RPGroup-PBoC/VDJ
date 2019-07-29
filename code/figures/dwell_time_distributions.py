#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import vdj.io 
import vdj.viz
vdj.viz.plotting_style()

#%%
# load the dwell time data
data = pd.read_csv('../../data/compiled_dwell_times.csv')
data = data[(data['salt']=='Mg') & (data['hmgb1']==80)]

# Load the fitting statistics. 
samps = pd.read_csv('../../data/expon_samples.csv')
samps = samps[samps['salt']=='Mg']
stats = pd.read_csv('../../data/expon_summary.csv')
stats = stats[stats['salt']=='Mg']

# select only the mutants of interest
endogenous = ['WT12rss', 'V19-93', 'V5-43'] 
heptamer = ['12HeptC3T', '12HeptA4T', '12HeptG7A']
spacer = ['12SpacC4T', '12SpacG11T', '12SpacA12C']
nonamer = ['12NonA1G', '12NonA3C', '12NonC8G']

rows = [endogenous, heptamer, spacer, nonamer]
colors = ['gray', 'tomato', 'dodgerblue', 'rebeccapurple']

DEADFILTER = 21 / 60
time = np.linspace(0, 60, 500) -  DEADFILTER


#%%
# Set up the figure canvas
fig, ax = plt.subplots(4, 3, figsize=(4.5, 5), sharex=True, sharey=True)

# Format the axes
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xscale('log')
    a.set_ylim([0, 1.1])
    a.set_xlim([DEADFILTER, 80])

# Add axis labels
for i in range(3):
    ax[-1,i].set_xlabel('dwell time [min]', fontsize=8)
for i in range(4):
    ax[i, 0].set_ylabel('cumulative\ndistribution', fontsize=8)

# Add row labels. 
fig.text(-0.05, 0.73, 'Endogenous', rotation='vertical', color='slategrey', 
         fontsize=9)
fig.text(-0.05, 0.55, 'Heptamer', rotation='vertical', color='tomato', 
         fontsize=9)
fig.text(-0.05, 0.37, 'Spacer', rotation='vertical', color='dodgerblue', 
         fontsize=9)
fig.text(-0.05, 0.16, 'Nonamer', rotation='vertical', color='rebeccapurple', 
         fontsize=9)
fig.text(-0.12, 0.87, '(A)', fontsize=9)
fig.text(-0.12, 0.67, '(B)', fontsize=9)
fig.text(-0.12, 0.47, '(C)', fontsize=9)
fig.text(-0.12, 0.27, '(D)', fontsize=9)



# Plot the observed dwell times. 
for i, row in enumerate(rows):
    for j, mut in enumerate(row):
        # Isolate the mutant
        mut_data = data[data['mutant']==mut]
        mut_stats = stats[stats['mutant']==mut]
        low = mut_stats[mut_stats['parameter']=='tau']['hpd_min'].values[0]
        high = mut_stats[mut_stats['parameter']=='tau']['hpd_max'].values[0]

        # Compute the empirical CDF
        x = list(np.sort(mut_data['dwell_time_min'] - DEADFILTER))
        y = list(np.arange(0, len(x), 1) / len(x))
        x.append(x[-1])
        y.append(1)
        # Plot the steps.
        ax[i, j].step(x, y, color='k', lw=1)

        # Compute the fit
        low= 1 - np.exp(-(time - DEADFILTER)/low)
        high = 1 - np.exp(-(time - DEADFILTER)/high)
        ax[i, j].fill_between(time, low, high, color=colors[i], alpha=0.5)

        # Add the titles
        if mut == 'WT12rss':
            title = 'V4-57-1'
        if 'Hept' in mut:
            title = f"Heptamer {mut.split('Hept')[-1]}"
        elif 'Spac' in mut:
            title = f'Spacer {mut.split("Spac")[-1]}'
        elif 'Non' in mut:
            title = f'Nonamer {mut.split("Non")[-1]}'
        else:
            title = mut
        ax[i, j].set_title(title, fontsize=8, y=0.95)

ax[0,0].set_title('V4-57-1', fontsize=8, y=0.95)
plt.savefig('../figures/FigX_nonexponential_distributions.pdf', 
            facecolor='white', bbox_inches='tight')

#%%
