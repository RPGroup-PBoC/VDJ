# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vdj.io
import vdj.viz
import scipy.stats as st
vdj.viz.plotting_style()

# Load the data and sampling information for single exponential
data = pd.read_csv('../../data/compiled_dwell_times.csv')
single_exp = pd.read_csv('../../data/expon_summary.csv')

#%%
def plot_ca_exp(data, single_exp_data, mutant, hmgb1=80, saveplot=False):
        single_ca = single_exp[(single_exp['mutant']==mutant) & (single_exp['hmgb1']==hmgb1)]
        mut = data[(data['mutant']==mutant) & (data['hmgb1']==hmgb1) & (data['salt']=='Ca')]
        
        ca_x, ca_y = np.sort(mut['dwell_time_min']), np.arange(0, len(mut), 1) / len(mut)
        time_range = np.linspace(0, 35, 500)
        t_min = 21/60

        cdf_min = 1 - np.exp(-(time_range - t_min)/single_ca['hpd_min'].values[0])
        cdf_max = 1 - np.exp(-(time_range - t_min)/single_ca['hpd_max'].values[0])
        
        plt.figure(dpi=100)
        if mutant=='WT12rss':
                plt.step(ca_x, ca_y, label='reference, 80 nM HMGB1', color='#753F98')
        else:
                plt.step(ca_x, ca_y, label='%s, 80 nM HMGB1' % mutant, color='#753F98')
        plt.fill_between(time_range, cdf_min, cdf_max, alpha=0.5, color='tomato',
                        label='single exponential')

        plt.ylim([0, 1.1])
        plt.title('Calcium')
        plt.legend(loc='lower right', fontsize=12)
        if saveplot:
                if mutant=='WT12rss':
                        plt.savefig('ref_ca_ecdf_posterior.pdf', facecolor='white')
                else:
                        plt.savefig('%s_ca_ecdf_posterior.pdf' % mutant, facecolor='white')
#%%
plot_ca_exp(data, single_exp, '12SpacG11T') 
#%%
def plot_ca_mg_ecdfs(data, mutant, hmgb1=80, saveplot=False):
        mut = data[(data['mutant']==mutant) & (data['hmgb1']==hmgb1)]
        mg = mut[mut['salt']=='Mg']
        ca = mut[mut['salt']=='Ca']

        ca_x, ca_y = np.sort(ca['dwell_time_min']), np.arange(0, len(ca), 1) / len(ca)
        mg_x, mg_y = np.sort(mg['dwell_time_min']), np.arange(0, len(mg), 1) / len(mg)

        plt.step(ca_x, ca_y, label='Calcium', color='#429447')
        plt.step(mg_x, mg_y, label='Magnesium', color='#753F98')

        plt.ylim([0, 1.1])
        plt.title('%s' % mutant)
        plt.legend(loc='lower right', fontsize=12)
        if saveplot:
                if mutant=='WT12rss':
                        plt.savefig('ref_ca_mg_ecdf.pdf', facecolor='white')
                else:
                        plt.savefig('%s_ca_mg_ecdf.pdf' % mutant, facecolor='white')
#%%
plot_ca_mg_ecdfs(data, '12HeptA4T', saveplot=True)
#%%
single_mg_spac = single_exp[(single_exp['mutant']=='12SpacG11T') & (single_exp['salt']=='Mg') & (single_exp['hmgb1']==80)]
single_ca_spac = single_exp[(single_exp['mutant']=='12SpacG11T') & (single_exp['salt']=='Ca') & (single_exp['hmgb1']==80)]

single_mg_ref = single_exp[(single_exp['mutant']=='WT12rss') & (single_exp['salt']=='Mg') & (single_exp['hmgb1']==80)]
single_ca_ref = single_exp[(single_exp['mutant']=='WT12rss') & (single_exp['salt']=='Ca') & (single_exp['hmgb1']==80)]

single_mg_hept = single_exp[(single_exp['mutant']=='12SpacG11T') & (single_exp['salt']=='Mg') & (single_exp['hmgb1']==80)]
single_ca_hept = single_exp[(single_exp['mutant']=='12SpacG11T') & (single_exp['salt']=='Ca') & (single_exp['hmgb1']==80)]

#%%
# Look only at one mutant
spac = data[(data['mutant']=='12SpacG11T') & (data['hmgb1']==80)]
mg_spac = spac[spac['salt']=='Mg']
ca_spac = spac[spac['salt']=='Ca']
#%%
ca_x_spac, ca_y_spac = np.sort(ca_spac['dwell_time_min']), np.arange(0, len(ca_spac), 1) / len(ca_spac)
time_range = np.linspace(0, 35, 500)
t_min = 21/60
# cdf = np.exp(-t_min/single_ca['median'].values[0]) - np.exp(-time_range/single_ca['median'].values[0])
cdf_min =  1 - np.exp(-(time_range - t_min)/single_ca_spac['hpd_min'].values[0])
cdf_max =  1 - np.exp(-(time_range - t_min)/single_ca_spac['hpd_max'].values[0])
plt.figure(dpi=100)
plt.step(ca_x_spac, ca_y_spac, label='12SpacG11T, 80 nM HMGB1', color='#753F98')
# plt.plot(time_range, cdf, color='tomato', label='single exponential')
plt.fill_between(time_range, cdf_min, cdf_max, alpha=0.5, color='tomato',
                label='single exponential')

plt.ylim([0, 1.1])
plt.title('Calcium')
plt.legend(loc='lower right',fontsize=12)
plt.savefig('12SpacG11T_ca_ecdf_posterior.pdf', facecolor='white')
#%%
ref = data[(data['mutant']=='WT12rss') & (data['hmgb1']==80)]
mg_ref = ref[ref['salt']=='Mg']
ca_ref = ref[ref['salt']=='Ca']

ca_x_ref, ca_y_ref = np.sort(ca_ref['dwell_time_min']), np.arange(0, len(ca_ref), 1) / len(ca_ref)
time_range = np.linspace(0, 35, 500)
t_min = 21/60
# cdf = np.exp(-t_min/single_ca['median'].values[0]) - np.exp(-time_range/single_ca['median'].values[0])
cdf_min =  1 - np.exp(-(time_range - t_min)/single_ca_ref['hpd_min'].values[0])
cdf_max =  1 - np.exp(-(time_range - t_min)/single_ca_ref['hpd_max'].values[0])
plt.figure(dpi=100)
plt.step(ca_x_ref, ca_y_ref, label='reference, 80 nM HMGB1', color='#753F98')
# plt.plot(time_range, cdf, color='tomato', label='single exponential')
plt.fill_between(time_range, cdf_min, cdf_max, alpha=0.5, color='tomato',
                label='single exponential')

plt.ylim([0, 1.1])
plt.title('Calcium')
plt.legend(loc='lower right',fontsize=12)
plt.savefig('ref12rss_ca_ecdf_posterior.pdf', facecolor='white')
#%%
hept = data[(data['mutant']=='12HeptA4T') & (data['hmgb1']==80)]
mg_hept = hept[hept['salt']=='Mg']
ca_hept = hept[hept['salt']=='Ca']

ca_x_hept, ca_y_hept = np.sort(ca_hept['dwell_time_min']), np.arange(0, len(ca_hept), 1) / len(ca_hept)
time_range = np.linspace(0, 35, 500)
t_min = 21/60
# cdf = np.exp(-t_min/single_ca['median'].values[0]) - np.exp(-time_range/single_ca['median'].values[0])
cdf_min =  1 - np.exp(-(time_range - t_min)/single_ca_hept['hpd_min'].values[0])
cdf_max =  1 - np.exp(-(time_range - t_min)/single_ca_hept['hpd_max'].values[0])
plt.figure(dpi=100)
plt.step(ca_x_hept, ca_y_hept, label='12HeptA4T, 80 nM HMGB1', color='#753F98')
# plt.plot(time_range, cdf, color='tomato', label='single exponential')
plt.fill_between(time_range, cdf_min, cdf_max, alpha=0.5, color='tomato',
                label='single exponential')

plt.ylim([0, 1.1])
plt.title('Calcium')
plt.legend(loc='lower right',fontsize=12)
plt.savefig('12HeptA4T_ca_ecdf_posterior.pdf', facecolor='white')

#%%
# Compare to magnesium
mg_x_spac, mg_y_spac = np.sort(mg_spac['dwell_time_min']), np.arange(0, len(mg_spac), 1) / len(mg_spac)
mg_x_ref, mg_y_ref = np.sort(mg_ref['dwell_time_min']), np.arange(0, len(mg_ref), 1) / len(mg_ref)
mg_x_hept, mg_y_hept = np.sort(mg_hept['dwell_time_min']), np.arange(0, len(mg_hept), 1) / len(mg_hept)

#%%


#%%
