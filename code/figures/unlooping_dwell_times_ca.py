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
        plt.xlabel('time (min)', fontsize=12)
        plt.ylabel('CDF', fontsize=12)
        plt.title('Calcium')
        plt.legend(loc='lower right', fontsize=12)
        if saveplot:
                if mutant=='WT12rss':
                        plt.savefig('ref_ca_ecdf_posterior.pdf', facecolor='white')
                else:
                        plt.savefig('%s_ca_ecdf_posterior.pdf' % mutant, facecolor='white')
#%%
plot_ca_exp(data, single_exp, '12HeptA4T') 
#%%
def plot_ca_mg_ecdfs(data, mutant, hmgb1=80, saveplot=False):
        mut = data[(data['mutant']==mutant) & (data['hmgb1']==hmgb1)]
        mg = mut[mut['salt']=='Mg']
        ca = mut[mut['salt']=='Ca']

        ca_x, ca_y = np.sort(ca['dwell_time_min']), np.arange(0, len(ca), 1) / len(ca)
        mg_x, mg_y = np.sort(mg['dwell_time_min']), np.arange(0, len(mg), 1) / len(mg)

        plt.step(ca_x, ca_y, label='Calcium', color='#753F98')
        plt.step(mg_x, mg_y, label='Magnesium', color='#429447')

        plt.ylim([0, 1.1])
        plt.xlabel('time (min)', fontsize=12)
        plt.ylabel('empirical CDF', fontsize=12)
        plt.title('%s' % mutant)
        plt.legend(loc='lower right', fontsize=12)
        if saveplot:
                if mutant=='WT12rss':
                        plt.savefig('ref_ca_mg_ecdf.pdf', facecolor='white')
                else:
                        plt.savefig('%s_ca_mg_ecdf.pdf' % mutant, facecolor='white')
#%%
plot_ca_mg_ecdfs(data, '12SpacG11T', saveplot=True)
#%%
ca_data = data[(data['hmgb1']==80) & (data['salt']=='Ca')]

ref = ca_data[ca_data['mutant']=='WT12rss']
hept = ca_data[ca_data['mutant']=='12HeptA4T']
spac = ca_data[ca_data['mutant']=='12SpacG11T']
#%%
ref_x, ref_y = np.sort(ref['dwell_time_min']), np.arange(0, len(ref), 1) / len(ref)
hept_x, hept_y = np.sort(hept['dwell_time_min']), np.arange(0, len(hept), 1) / len(hept)
spac_x, spac_y = np.sort(spac['dwell_time_min']), np.arange(0, len(spac), 1) / len(spac)

plt.figure(dpi=100)
plt.step(ref_x, ref_y, label='reference', color='#753F98')
plt.step(hept_x, hept_y, label='12HeptA4T', color='#429447')
plt.step(spac_x, spac_y, label='12SpacG11T', color='#D43124')

plt.ylim([0, 1.1])
plt.title('Calcium, 80 nM HMGB1')
plt.xlabel('time (min)', fontsize=12)
plt.ylabel('empirical CDF', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.savefig('ca_ecdf_comparisons.pdf', facecolor='white')
#%%
