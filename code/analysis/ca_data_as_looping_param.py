# -*- coding: utf-8 -*-
#%% 
import numpy as np
import pandas as pd
import pystan
import matplotlib.pyplot as plt
import scipy.stats
import vdj.io
import vdj.viz
import glob
vdj.viz.plotting_style()



# %%
# Load th data
files = glob.glob('../../data/12SpacG11T*.mat')

dwell, floop, cuts = [], [], []
for f in files:
    tpm = vdj.io.ProcessTPM(f)
    f_loop, dwell_t, fates = tpm.extract_data()

    if '_ca_' in f.lower():
        salt = 'Ca'
    else:
        salt = 'Mg'
    
    f_loop['salt'] = salt
    fates['salt'] = salt
    dwell_t['salt'] = salt
    dwell.append(dwell_t)
    floop.append(f_loop)
    cuts.append(fates) 

# %%
dwell = pd.concat(dwell)
floop = pd.concat(floop)
cuts = pd.concat(cuts)

#%%

model_code = """

data { 
    int<lower=0> N_mg; 
    int<lower=0> N_ca;
    vector<lower=21>[N_mg] dwell_mg;
    vector<lower=21>[N_ca] dwell_ca;

    int<lower=0> n_cuts;
    int<lower=0> n_beads;
    real<lower=0, upper=1> f_loop;
}

parameters { 
    real<lower=0> k_loop;
    real<lower=0> k_unloop;
    real<lower=0> r_cut;
}

model {
    k_loop ~ normal(0, 1);
    k_unloop ~ normal(0, 1);
    r_cut ~ normal(0, 1);

    dwell_ca ~ exponential(k_unloop);
    dwell_mg ~ exponential(k_unloop + r_cut);

    n_cuts ~ binomial(n_beads, r_cut / (r_cut + k_unloop));

}

generated quantities {
    real k_loop = (k_unloop + r_cut) / ((1 / f_loop) - 1);
    real p_loop = k_loop / (k_unloop + k_loop);
    real p_cut = r_cut / (r_cut + k_unloop);
    real f_loop = k_loop / (k_unloop + k_loop + r_cut);

}

"""
model = pystan.StanModel(model_code=model_code)

# %%

# Generate the data dict
dwell_ca = dwell[dwell['salt']=='Ca']
dwell_mg = dwell[dwell['salt']=='Mg']
fates_mg = cuts[cuts['salt']=='Mg']
floop_ca = floop[floop['salt']=='Ca']
floop_mg = floop[floop['salt']=='Mg']

# %%
data_dict = {'N_mg':len(dwell_mg), 'N_ca':len(dwell_ca),
             'dwell_mg':dwell_mg['dwell_time_s'],
             'dwell_ca':dwell_ca['dwell_time_s'],
             'n_cuts':fates_mg['n_cuts'].sum().astype(int),
             'n_beads':fates_mg['n_beads'].sum().astype(int),
             'n_looped_mg':floop_mg['looped_frames'].sum().astype(int),
             'n_frames_mg':floop_mg['total_frames'].sum().astype(int),
             'n_looped_ca':floop_ca['looped_frames'].sum().astype(int),
             'n_frames_ca':floop_ca['total_frames'].sum().astype(int),
             'f_loop':floop_mg['looped_frames'].sum() / floop_mg['total_frames']..sum()}

samples = model.sampling(data_dict)
samples_df = samples.to_dataframe()

# %%
samples
#%%

# Look at the dwell_times
fig, ax = plt.subplots(1, 2, figsize=(6, 3), dpi=200, sharex=True, sharey=True)
ax[0].hist(dwell_ca['dwell_time_s'], bins=100, density=True, edgecolor='dodgerblue')
ax[1].hist(dwell_mg['dwell_time_s'], bins=100, density=True,
edgecolor='tomato')

# Plot the exponential fits.
dwell_range = np.logspace(0, 4, 500) + 21
ca_exp = scipy.stats.expon(loc=0, scale=(1 / (samples_df['k_unloop'].mean()))).pdf(dwell_range)
mg_exp = scipy.stats.expon(loc=0, scale=(1 / (samples_df['r_cut'].mean() + samples_df['k_unloop'].mean()))).pdf(dwell_range)
ax[1].plot(dwell_range, ca_exp, 'k-')
ax[0].plot(dwell_range, mg_exp, 'k-')

for a in ax:
    a.set_xlabel('dwell time [s]')
    a.set_ylabel('frequency')

ax[0].set_title('with magnesium')
ax[1].set_title('with calcium')


plt.tight_layout()
#%% [markdown]

#%%

# Plot the ecdfs
ca_x  = np.sort(dwell_ca['dwell_time_s'].values) 
ca_y = np.linspace(0, 1, len(ca_x))
mg_x  = np.sort(dwell_mg['dwell_time_s'].values) 
mg_y = np.linspace(0, 1, len(mg_x))


# %%
fig, ax = plt.subplots(1,1, dpi=200)
ax.step(ca_x, ca_y, label='Ca$^{2+}$', color='tomato')
# ax.step(mg_x, mg_y, label='Mg$^{2+}$', color='dodgerblue')

ca_cdf = 1 - np.exp(-samples_df.k_unloop.median() * (dwell_range))
mg_cdf = 1 - np.exp(-(samples_df.k_unloop.median() + samples_df.r_cut.median())* (dwell_range))
ax.step(dwell_range, ca_cdf,  color='tomato', lw=1, linestyle=':', label='Ca fit')
# ax.step(dwell_range, mg_cdf,  color='dodgerblue', lw=1, linestyle=':', label='Mg fit')
ax.set_xlabel('dwell time [s]')
ax.set_ylabel('ECDF')
ax.legend()

ax.set_xscale('log')




##%%

dwell_ca


#%%
plt.hist(samples_df['p_cut'])

#%%

# Look only at the calcium case. 

model_code = """
data { 
int N;
vector[N] dwell;
}

parameters {
    real k;
}

model { 
    k ~ normal(0, 0.1);
    dwell ~ exponential(k);
}
"""

ca_model = pystan.StanModel(model_code=model_code)

#%%
data_dict = {'N':len(dwell_ca), 'dwell':dwell_ca.dwell_time_s - 21}
ca_samples = ca_model.sampling(data=data_dict)

#%%
ca_df = ca_samples.to_dataframe()
cred_region = np.zeros((2, len(dwell_range)))
for i, t in enumerate(dwell_range):
    ca_cdf = 1 - np.exp(-ca_df['k'] * (t-21))
    cred_region[:, i] = vdj.stats.compute_hpd(ca_cdf, 0.99)

plt.figure(dpi=200)
plt.step(ca_x-21, ca_y, 'k', label = 'Ca$^{2+}$ data')
plt.fill_between(dwell_range, cred_region[0, :], cred_region[1, :],
color='tomato', alpha=0.5, label='$1 - e^{-k_{ul}(t - t_{min})}$')
plt.xscale('log')
plt.legend(fontsize=12)
plt.xlabel('dwell time [s]')
plt.ylabel('cumulative distribution')
plt.savefig('/Users/gchure/ca_data_dwell.pdf', bbox_inches='tight', transparent=False)
#%%
model_code = """
data {
    int N;
    vector[N] dwell;
    int N_l;
    int n_c; 
}

parameters {
    real k_ul;
    real r_cut;
}

model  {
    k_ul ~ normal(0, 0.1);
    r_cut ~ normal(0, 0.1);
    dwell ~ exponential(k_ul + r_cut);
    n_c ~ binomial(N_l, r_cut / (k_ul + r_cut));
}
"""
mg_model = pystan.StanModel(model_code=model_code)

#%%
data_dict = {'N': len(dwell_mg), 'dwell':dwell_mg.dwell_time_s - 21, 
            'N_l':fates_mg['n_beads'].sum(), 'n_c':fates_mg['n_cuts'].sum()}

mg_samples = mg_model.sampling(data_dict)

#%%
mg_df = mg_samples.to_dataframe()
mg_cdf = 1 - np.exp(-(mg_df.k_ul.mean() + mg_df.r_cut.mean()) * (dwell_range -21))


plt.figure(dpi=200)
plt.step(mg_x-21, mg_y, label = 'Mg$^{2+}$ data', color='tomato')
plt.step(dwell_range, mg_cdf,  'k-', label='$1 - e^{-(k_{ul} + r_{cut}) * (t - t_{min})}$')
plt.xscale('log')
plt.legend(fontsize=12)
plt.xlabel('dwell time [s]')
plt.ylabel('cumulative distribution')
plt.savefig('/Users/gchure/ca_data_dwell.pdf')


#%%
model_code = """
functions {
    real dwell_dist_lpdf(vector dwell, vector rates) {
        return sum(log(rates[1]*exp(-rates[1] * dwell) + rates[2]*exp(-rates[2] * dwell)));
    }
}

data {
    int N;
    vector[N] dwell;
}

parameters {
   positive_ordered[2] k; 
}

model {
    k ~ normal(0, 1);
    dwell ~ dwell_dist(k);
}
"""
tworate_model = pystan.StanModel(model_code=model_code)
#%%
data_dict = {'N':len(dwell_ca), 'dwell':dwell_ca.dwell_time_s}   
tworate_samples = tworate_model.sampling(data_dict);
tworate_df = tworate_samples.to_dataframe()

#%%
mg_df = mg_samples.to_dataframe()
k1 = tworate_df['k[1]'].mean()
k2 = tworate_df['k[2]'].mean()
pdf = k1 * np.exp(-k1 * (dwell_range - 21)) + k2 * np.exp(-k2 * (dwell_range - 21))
plt.plot(dwell_range , pdf, color='tomato', lw=2)


plt.hist(dwell_ca.dwell_time_s - 21, bins=50, density=True, color='dodgerblue', lw=0)
# plt.xscale('log')

plt.figure()


#%%
