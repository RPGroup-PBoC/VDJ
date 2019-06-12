# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import vdj.io
import vdj.stats
import vdj.viz
import seaborn as sns
colors = sns.color_palette('deep')
vdj.viz.plotting_style()

# Load the data and sampling info
data = pd.read_csv('../../data/compiled_dwell_times.csv')
sum_samps = pd.read_csv('../../data/sum_expon_samples.csv')
mono_samps = pd.read_csv('../../data/expon_samples.csv')
mono_stats = pd.read_csv('../../data/expon_summary.csv')
# %% Isolate data and samples to only the calcium data.
ca_data = data[data['salt'] == 'Ca']
ca_sum = sum_samps[sum_samps['salt']== 'Ca']
ca_mono = mono_samps[mono_samps['salt']=='Ca']
ca_mono_stats = mono_stats[mono_stats['salt']=='Ca']
# %%
# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=100)

# Plot the ecdf of the dwell times
ax.set_xlabel('paired complex dwell time [min]', fontsize=12)
ax.set_ylabel('frequency', fontsize=12)
hist, bins = np.histogram(ca_data['dwell_time_min'], bins=100)
hist = hist / np.sum(hist)
ax.step(bins[:-1], hist, label='data', color=colors[0], lw=0.75)
ax.fill_between(bins[:-1], hist, label='data', color=colors[0], alpha=0.5)
# ax.step(x, y, '.', color=colors[0], label='data', ms=4, markerfacecolor='w', 
        # mew=0.75)
# ax.set_xscale('log')

# Set up the time range and evaluate the credible region
time_range = np.linspace(0, 45, 800)
mono_pdf = scipy.stats.expon(loc=0, scale=(1/ca_mono_stats['median'])).pdf(time_range) 
mono_pdf = mono_pdf / np.sum(mono_pdf)
ax.plot(time_range, mono_pdf, color='tomato')
# ax.set_xscale('log')
# mono_cred_region =  n.zeros((2, len(time_range)))
# sum_cred_region = np.zeros((2, len(time_range)))
# for i, t in enumerate(time_range):
#     pdf = (1/ca_mono['tau']) * np.exp(-t/ca_mono['tau'])
#     mono_cred_region[:, i] = vdj.stats.compute_hpd(cdf, 0.95)
#     pdf = (ca_sum['theta'] / ca_sum['tau[1]'])*np.exp(-t/ca_sum['tau[1]']) + ((1 - ca_sum['theta'])/ca_sum['tau[2]']) * np.exp(-t/ca_sum['tau[2]'])
#     sum_cred_region[:, i] = vdj.stats.compute_hpd(cdf, 0.95)

# ax.fill_between((time_range) * 60, mono_cred_region[0, :] / np.sum(mono_cred_region[0, :]), mono_cred_region[1, :] / np.sum(mono_cred_region[1, :]), color=colors[1], alpha=0.7, label='exponential')
# ax.fill_between((time_range) * 60, sum_cred_region[0, :], sum_cred_region[1, :], color=colors[2], alpha=0.7, label='sum of exponentials')

# ax.set_title('Ca$^{2+}$ dwell times for 12SpacG11T', fontsize=12)
# ax.set_xscale('log')
# ax.set_ylim([-0.01, 1.01])
# ax.legend()
# plt.savefig('./ca_dwell_time_models.pdf', bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
#%%

# Now look at the case for Mg 2+

# %% Isolate data and samples to only the calcium data.
mg_data = data[(data['salt'] == 'Mg') & (data['mutant']=='12SpacG11T')]
mg_sum = sum_samps[(sum_samps['salt']== 'Mg') &\
         (sum_samps['mutant']=='12SpacG11T')]
mg_mono = mono_samps[(mono_samps['salt']=='Mg') &\
             (mono_samps['mutant']=='12SpacG11T')]

#%%
# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=100)

# Plot the ecdf of the dwell times
x = np.sort(mg_data['dwell_time_min'])  * 60
y = np.arange(0, len(mg_data)) / len(mg_data)
ax.set_xlabel('paired complex dwell time [sec]', fontsize=12)
ax.set_ylabel('cumulative distribution', fontsize=12)
ax.step(x, y, label='__nolegend__', color=colors[0], lw=0.75)
ax.step(x, y, '.', color=colors[0], label='data', ms=4, markerfacecolor='w', 
        mew=0.75)
# ax.set_xscale('log')

# Set up the time range and evaluate the credible region
time_range = np.linspace(0, 45, 800)
mono_cred_region =  np.zeros((2, len(time_range)))
sum_cred_region = np.zeros((2, len(time_range)))
for i, t in enumerate(time_range):
    cdf = 1 - np.exp(-t/mg_mono['tau'])
    mono_cred_region[:, i] = vdj.stats.compute_hpd(cdf, 0.95)
    cdf = -mg_sum['theta']*np.exp(-t/mg_sum['tau[1]']) + (mg_sum['theta'] - 1) * np.exp(-t/mg_sum['tau[2]']) + 1
    sum_cred_region[:, i] = vdj.stats.compute_hpd(cdf, 0.95)

ax.fill_between((time_range) * 60, mono_cred_region[0, :], mono_cred_region[1, :], color=colors[1], alpha=0.7, label='exponential')
ax.fill_between((time_range) * 60, sum_cred_region[0, :], sum_cred_region[1, :], color=colors[2], alpha=0.7, label='sum of exponentials')

ax.set_title('Mg$^{2+}$ dwell times for 12SpacG11T', fontsize=12)
ax.set_xscale('log')
ax.set_ylim([-0.01, 1.01])
ax.legend()
plt.savefig('./mg_dwell_time_models.pdf', bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')

#%%


#%%
