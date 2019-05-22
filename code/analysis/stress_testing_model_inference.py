#%% [markdown]
# # Stress Testing the Inferrential Model
# ---

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import vdj.viz
import vdj.bayes
import vdj.stats
import imp
imp.reload(vdj.viz)
vdj.viz.plotting_style()

#%% [markdown]
# In this notebook, we play around with the inferential scheme for the model
# parameters we are interested in to determine if they are a) identifiable and
# b) to determine what the diminishing returns are on doing replicates.
#
# For right, now I will not go through the details of the model and will simply
# implement it here.

# ## Assigning prior distributions
# Rather than doing the complete principled workflow, let's just choose "true"
# values for the parameters and distribute some data about them with different
# numbers of observations.
#%%
# Define the ground-truth parameters
# Define the "true" model parameters
r_cut_mu = 1 
r_cut_sigma = 0.1 
k_loop_mu = 5
k_loop_sigma = 0.1
k_unloop_mu = 1 
k_unloop_sigma = 0.1
r_mu = r_cut_mu + k_unloop_mu
sigma_truth = 0.001

# Define the number of replicates and generate a storage vector for dfs
n_rep = 100 
dfs = []
for n in range(n_rep): 
    # Determine the number of beads for that day
    n_beads = int(np.random.normal(12, 1.5))

    # Define the data sets
    r_cut = np.random.normal(r_cut_mu, r_cut_sigma);
    k_unloop = np.random.normal(k_unloop_mu, k_unloop_sigma); 
    k_loop = np.random.normal(k_loop_mu, k_loop_sigma);
    dwell_time = np.random.exponential(r_cut + k_unloop, n_beads)
    cut = [np.random.binomial(1, r_cut/(r_cut + k_unloop)) for _ in range(n_beads)]
    f_loop = np.random.normal((k_loop / (r_cut + k_unloop + k_loop)), sigma_truth)
    _df = pd.DataFrame([])
    _df['dwell_time'] = dwell_time
    _df['cut'] = cut
    _df['f_loop'] = f_loop
    _df['rep'] = n+1
    dfs.append(_df)

# Assemble into a single dataframe
df = pd.concat(dfs)

#%%[markdown]

# ## Model 1: Pooled data
# The first and simplest approach to take is to assume that there is no
# day-to-day variation in the rates or sample preparation, or at least the error
# introduced by the variation is far below the noise level of the system. 

#%%
# Load the stan model
model = vdj.bayes.StanModel('../stan/pooled_model.stan', force_compile=True)
 
#%%
# Define the data dictionary
data_dict = {'N':len(df), 'cut':df['cut'], 'f_looped':np.mean(df['f_loop'].values),
            'dwell_time':df['dwell_time']}
samples, samples_df = model.sample(data_dict, iter=5000)

#%%
samples

#%%
x, y = np.sort(df['dwell_time']), np.linspace(0, 1, len(df))
plt.plot(x, y, 'k-')
x_range = np.linspace(0, 1.5, 500)
exp_cdf = scipy.stats.expon(1/2.55).cdf(x_range)
plt.plot(x_range, exp_cdf, '-')

#%%
samples_df.head()

#%%
# ## Model II: Hierarchical
# The proper way to do this is as a hierarchical model in which the day-to-day
# variation is taken into account.
#%%
model = vdj.bayes.StanModel('../stan/hierarchical_model.stan', force_compile=True)

#%%
# Assemble the data dictionary. 
data_dict = {'J':df['rep'].max(), 'N':len(df), 'idx':df['rep'].values,
             'f_looped':df['f_loop'].unique(), 'cut':df['cut'], 'dwell_time':df['dwell_time']}
samples, samples_df = model.sample(data_dict, iter=6000, control=dict(adapt_delta=0.99))

#%%
samples


#%%[markdown]
# ## Trying with some data

# Let's try wth some actual data. I'm not convinced that I'm generating the data
# with believable values. Guess I don't have the domain expertise!
import scipy.io
mat = scipy.io.loadmat('../../data/analysis_280_analyzed.mat')
#%%
mat.keys()


#%%
