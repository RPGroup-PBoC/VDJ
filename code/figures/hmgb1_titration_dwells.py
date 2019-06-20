# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import vdj.io
import vdj.viz
from tqdm import tqdm
import scipy.stats as st
vdj.viz.plotting_style()

# Load dwell time data
dwell_data = pd.read_csv('../../data/compiled_dwell_times.csv')

# Pull out 12SpacG11T data
hmgb1_titration_data = dwell_data[(dwell_data['mutant']=='12SpacG11T') & (dwell_data['salt']=='Ca')]
#%%
# Plot ECDFs for each HMGB1 concentration
dwell_40 = hmgb1_titration_data[hmgb1_titration_data['hmgb1']==40]
dwell_80 = hmgb1_titration_data[hmgb1_titration_data['hmgb1']==80]
dwell_160 = hmgb1_titration_data[hmgb1_titration_data['hmgb1']==160]

x_40, y_40 = np.sort(dwell_40['dwell_time_min']), np.arange(1, len(dwell_40)+1, 1) / len(dwell_40)
x_80, y_80 = np.sort(dwell_80['dwell_time_min']), np.arange(1, len(dwell_80)+1, 1) / len(dwell_80)
x_160, y_160 = np.sort(dwell_160['dwell_time_min']), np.arange(1, len(dwell_160)+1, 1) / len(dwell_160)
# %%
# Bootstrapping 40 nM HMGB1 dwell times
bootstrap_num = 100000
sample_40 = np.random.choice(dwell_40['dwell_time_min'], 
                                size=(bootstrap_num,len(dwell_40)), 
                                replace=True)
sample_160 = np.random.choice(dwell_160['dwell_time_min'],
                                size=(bootstrap_num, len(dwell_160)),
                                replace=True)
x_40_sample, y_40_sample = np.sort(sample_40,axis=1), np.arange(1, len(dwell_40)+1, 1) / len(dwell_40)
x_160_sample, y_160_sample = np.sort(sample_160, axis=1), np.arange(1, len(dwell_160)+1, 1) / len(dwell_160)

#Compute 25th to 75th percentiles of x_sample values for each CDF value
bootstrap_40_25 = np.percentile(x_40_sample, 25, axis=0)
bootstrap_40_75 = np.percentile(x_40_sample, 75, axis=0)

bootstrap_160_25 = np.percentile(x_160_sample, 25, axis=0)
bootstrap_160_75 = np.percentile(x_160_sample, 75, axis=0)
#%%
plt.figure(dpi=100)
plt.fill_betweenx(y_40_sample, bootstrap_40_25, bootstrap_40_75, 
                step='pre', facecolor='#429447',
                alpha=0.65, label='40 nM bootstrap, %i loops' %len(dwell_40))
plt.fill_betweenx(y_160_sample, bootstrap_160_25, bootstrap_160_75,
                step='pre', facecolor='#D43124',
                alpha=0.65, label='160 nM bootstrap, %i loops' %len(dwell_160))
plt.step(x_80, y_80, label='80 nM, %i loops' %len(dwell_80), color='#753F98')
plt.xlabel('time (min)')
plt.ylabel('bootstrapped ECDFs')
plt.legend(loc='lower right')
plt.title('HMGB1 Titration (Bootstrapping)')
plt.savefig('hmgb1_titration_bootstrap.pdf', facecolor='w')
#%%
plt.figure(dpi=100)
plt.step(x_40, y_40, label='40 nM, %i loops' %len(dwell_40), color='#429447')
plt.step(x_80, y_80, label='80 nM, %i loops' %len(dwell_80), color='#753F98')
plt.step(x_160, y_160, label='160 nM, %i loops' %len(dwell_160), color='#D43124')
plt.xlabel('time (min)')
plt.ylabel('empirical CDF')
plt.title('HMGB1 Titration (Calcium)')
plt.legend()
plt.savefig('hmgb1_titration_raw_dwells.pdf', facecolor='w')
#%%
