# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.io
import vdj.viz
import scipy.stats as st
vdj.viz.plotting_style()

# Load the data and sampling information
data = pd.read_csv('../../data/compiled_dwell_times.csv')
single_exp = pd.read_csv('../../data/expon_summary.csv')
sum_exp = pd.read_csv('../../data/sum_expon_summary.csv')
single_mg = single_exp[(single_exp['mutant']=='12SpacG11T') & (single_exp['salt']=='Mg')]
single_ca = single_exp[(single_exp['mutant']=='12SpacG11T') & (single_exp['salt']=='Ca')]
sum_mg = sum_exp[(sum_exp['mutant']=='12SpacG11T') & (sum_exp['salt']=='Mg')]
sum_ca = sum_exp[(sum_exp['mutant']=='12SpacG11T') & (sum_exp['salt']=='Ca')]

#%%
# Look only at one mutant
mut = data[data['mutant']=='12SpacG11T']
mg = mut[mut['salt']=='Mg']
ca = mut[mut['salt']=='Ca']
#%%
# Look at the distributions
tau1 = sum_mg[(sum_mg['parameter']=='tau1')]['median'].values[0]
tau2 = sum_mg[(sum_mg['parameter']=='tau2')]['median'].values[0]
theta = sum_mg[(sum_mg['parameter']=='theta')]['median'].values[0]

time_range = np.linspace(0, 50, 500)
dist = (1/single_mg['median'].values[0]) * np.exp(-time_range / single_mg['median'].values[0])
sum_dist = (theta/tau1) * np.exp(-time_range/tau1) + ((1 - theta) / tau2) * np.exp(-time_range/tau2)
fig, ax = plt.subplots(1,1, dpi=100)
ax.hist(mg['dwell_time_min'], bins=50, color='grey', alpha=0.5, density=True)
ax.plot(time_range, dist, color='tomato', label='single exponential')
ax.plot(time_range, sum_dist, color='dodgerblue', label='sum of exponential')

# %%
# Look at the distributions
tau1 = sum_ca[(sum_ca['parameter']=='tau1')]['median'].values[0]
tau2 = sum_ca[(sum_ca['parameter']=='tau2')]['median'].values[0]
theta = sum_ca[(sum_ca['parameter']=='theta')]['median'].values[0]

time_range = np.linspace(0, np.max(ca['dwell_time_min']), 500)
dist = (1/single_ca['median'].values[0]) * np.exp(-time_range / single_ca['median'].values[0])
sum_dist = (theta/tau1) * np.exp(-time_range/tau1) + ((1 - theta) / tau2) * np.exp(-time_range/tau2)
# sum_dist = sum_dist / np.sum(sum_dist)
# dist = dist / np.sum(dist)
hist, bins = np.histogram(ca['dwell_time_min'], bins=75)
hist = hist / np.sum(hist)
fig, ax = plt.subplots(1,1, dpi=100)
ax.step(bins[:-1], hist, color='grey')
ax.fill_between(bins[:-1], hist, color='grey', alpha=0.5, step='pre')
# ax.hist(ca['dwell_time_min'], bins=50, color='grey', alpha=0.5, density=True)
ax.plot(time_range, dist, color='tomato', label='single exponential')
ax.plot(time_range, sum_dist, color='dodgerblue', label='sum of exponential')



#%%
ca_x, ca_y = np.sort(ca['dwell_time_min']), np.arange(0, len(ca), 1) / len(ca)
time_range = np.linspace(0, 35, 500)
t_min = 21/60
# cdf = np.exp(-t_min/single_ca['median'].values[0]) - np.exp(-time_range/single_ca['median'].values[0])
cdf_min =  1 - np.exp(-(time_range - t_min)/single_ca['hpd_min'].values[0])
cdf_max =  1 - np.exp(-(time_range - t_min)/single_ca['hpd_max'].values[0])
plt.figure(dpi=100)
plt.step(ca_x, ca_y, label='data, 80 nM HMGB1', color='#753F98')
# plt.plot(time_range, cdf, color='tomato', label='single exponential')
plt.fill_between(time_range, cdf_min, cdf_max, alpha=0.5, color='tomato',
                label='single exponential')

plt.ylim([0, 1.1])
plt.title('Calcium')
plt.legend(loc='lower right',fontsize=12)
plt.savefig('12SpacG11T_ca_ecdf_posterior.pdf',
            facecolor='white')
#%%
time_range =np.linspace(0, 50, 500)
mg_x, mg_y = np.sort(mg['dwell_time_min']), np.arange(0, len(mg), 1) / len(mg)
tau1 = sum_mg[(sum_mg['parameter']=='tau1')]['median'].values[0]
tau2 = sum_mg[(sum_mg['parameter']=='tau2')]['median'].values[0]
theta = sum_mg[(sum_mg['parameter']=='theta')]['median'].values[0]


t_min = 21/60
# cdf = np.exp(-t_min/single_ca['median'].values[0]) - np.exp(-time_range/single_ca['median'].values[0])
cdf_min =  1 - np.exp(-(time_range - t_min)/single_mg['hpd_min'].values[0])
cdf_max =  1 - np.exp(-(time_range - t_min)/single_mg['hpd_max'].values[0])

numer = theta * np.exp(-time_range/tau1) + (1 - theta) * np.exp(-time_range / tau2)
denom = theta * np.exp(-t_min/tau1) + (1 - theta) * np.exp(-t_min / tau2)

dcdf = 1 - (numer / denom) 
fig = plt.figure(dpi=100)
plt.step(mg_x, mg_y, label='ECDF', color='#753F98')
# plt.plot(time_range, cdf, color='tomato', label='single exponential')
plt.fill_between(time_range, cdf_min, cdf_max, color='tomato', 
                label=r'$\frac{1}{\tau}e^{-(t-t_{min})/\tau}$ pdf', alpha=0.5)
plt.plot(time_range, dcdf, 
        label=r'$\approx\frac{p}{\tau_1}e^{-t/\tau_1} + \frac{(1-p)}{\tau_2}e^{-t/\tau_2}$ pdf', 
        color='dodgerblue')
plt.ylim([0, 1.1])
plt.title('Magnesium')
plt.legend(loc='lower right',fontsize=12)
plt.savefig('12SpacG11T_mg_ecdf_single_double_exp_posterior.pdf', 
            facecolor='white')

#%%
