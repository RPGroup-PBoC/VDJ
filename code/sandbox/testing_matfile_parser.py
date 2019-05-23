# -*- coding: utf-8 -*-
# %%[markdown]
# # Developing a TPM class
#%%
import numpy as np
import pandas as pd
import scipy.io
import vdj.io
import vdj.viz
import vdj.bayes
import matplotlib.pyplot as plt
from functools import reduce
import imp
imp.reload(vdj.bayes)
vdj.viz.plotting_style()

# In this notebook, we we test the matfile input-output reading function to make
# sure that everything is beging computed correctly


#%%
_mat = scipy.io.loadmat('../../data/analysis_280_analyzed.mat',
matlab_compatible=True)
_mat.keys()
# %%
# Definitely need to define a class to do all this stuff. 
class ProcessTPM(object):
    def __init__(self, fname=None, matfile=None, framerate=30, 
                 stan_model=None):
        if (fname == None) & (matfile != None):
            self.file = matfile
        elif (fname != None) & (matfile == None):
            self.file = scipy.io.loadmat(fname)
        else:
            raise RuntimeError("Either `fname` or `matfile` must be specified, but not both!")

        if framerate <= 0:
            raise ValueError('framerate must be greater than zero!')

        # Define variables
        self.framerate = framerate # in ms/frame
        self.stan_model = stan_model

        # Initialize the extracted quantities
        self.dwell = False
        self.f_looped = False
        self.fates = False

        # Set state variable for if data is processed
        self.extract = False

    # ##########################################################################
    # DATA EXTRACTION
    # ##########################################################################
    def dwell_time(self):
        """
        Extracts the dwell times for each replicate.
        """
        mat = self.file
        dfs = [] # empty list to append coming data frames
        iter = 1
        for i, rep in enumerate(mat['loops'][0]):
            # Find only the times where loopstate was observed
            dwell = rep[rep > 0]
            # If there were no looped states, increment the counter and move on
            if len(dwell) == 0:
                dwell = np.array([0])
                iter += 1
            else:
                _df = pd.DataFrame([])
                _df['dwell_time_ms'] = dwell 
                _df['replicate'] = iter # indexing replicates by 1 
                iter +=1
                dfs.append(_df)
        df = pd.concat(dfs).reset_index()
 
        # Make the appropriate entries integers
        df['replicate'] = df['replicate'].values.astype(int)
        self.dwell = df
        return df
 
    def fraction_looped(self):
        """
        Computes the fraction of time the beads were in the looped state for
        each replicate. This is defined as 
                    Σ τ_looped
        f_looped = ------------ 
                    Σ τ_total
        """
        mat = self.file
        df = pd.DataFrame([])
        for i, rep in enumerate(mat['statetrace_comp'][0]):
            # Determine the total length of the experiment and the number of
            # beads observed
            n_beads, total_time = np.shape(rep)
 
            # Identify the looped and unlooped states. Key for these idx are 
            # 1 = Bead stuck to glass or cut
            # 3 = Bead in unlooped state
            # 2 = Bead in looped state
            total_time = len(rep[rep != 1])
            f_looped = np.sum(rep == 2) / total_time
            df = df.append({'fraction_looped':f_looped,
                            'total_frames': total_time,
                            'looped_frames': np.sum(rep==2),
                            'n_beads': n_beads,
                            'replicate': i + 1}, ignore_index=True)
 
        # Make the appropriate entries integers
        df['n_beads'] = df['n_beads'].values.astype(int)
        df['replicate'] = df['replicate'].values.astype(int)
        df['total_frames'] = df['total_frames'].values.astype(int)
        df['looped_frames'] = df['looped_frames'].values.astype(int)
        self.f_looped = df.reset_index() 
        return df
 
    def cut_beads(self, bead_idx=False):
        """ 
        Generates a DataFrame of the replicates and bead id for each cutting
        event.
        
        Parameters
        ----------
        bead_idx: bool
            If True, the returned dataframe has the bead identity and dwell
            time. If False, a short dataframe with n_cuts and n_beads for each
            replicate are returned.
        """
        mat = self.file
        df = pd.DataFrame({})
        for i, rep in enumerate(mat['pccut'][0]):
            # Generate a dictionary of all dwell times and bead ids for each
            # bead in the replicate 
            bead_dict = {}
            n_beads = np.sum(mat['ontime_comp'][0][i] > 0)
            for k, bead in enumerate(mat['ontime_comp'][0][i]):
                if bead.sum() > 0:
                    dwells = bead[bead > 0]
                    for d in dwells:
                        bead_dict[d] = k + 1
 
            # Trim the replicate to remove spurious entry of zero
            rep = rep[0][1:]
 
            # Get the bead ID and add it to the dataframe
            if  bead_idx == True:
                for dwell in rep:
                    df = df.append({'cut':1, 
                                'bead_id': bead_dict[dwell],
                                'replicate': i + 1}, ignore_index=True)
            else:
                df = df.append({'n_beads':n_beads, 'n_cuts':len(rep), 
                                'replicate': i + 1}, ignore_index=True)
                
        # Make the appropriate entries integers
        if bead_idx == True:
            df['bead_id'] = df['bead_id'].values.astype(int)
            self.bead_idx = True
        else:
            df['n_beads'] = df['n_beads'].values.astype(int)
            df['n_cuts'] = df['n_cuts'].values.astype(int)
            self.bead_idx = False

        df['replicate'] = df['replicate'].values.astype(int)
        self.fates = df.reset_index()
        return df
 
    def extract_data(self, bead_idx=False):
        """
        Extracts the dwell times, fraction looped, and cutting identity and
        returns a list of the three data frames
        """
        if type(self.f_looped) == bool:
            f_looped = self.fraction_looped()
        else:
            f_looped = self.f_looped
        if type(self.dwell) == bool:
            dwell = self.dwell_time()
        else:
            dwell = self.dwell
        if type(self.fates) == bool:
            fates = self.cut_beads(bead_idx=bead_idx)
        else:
            fates = self.fates
        self.extract = True
        return [f_looped, dwell, fates]
 

    # ##########################################################################
    # INFERENCE 
    # ##########################################################################
    def run_inference(self, stan_model=None, force_compile=False, iter=2000, sampler_kwargs={}):
        """
        Executes the model inference using provided stan model. 

        Parameters
        ----------
        stan_model : str
            Path to the stan model. Note that this function *does not* take the
            model code as a string. the Stan model must be saved as a separate
            file!
        
        Returns [as list]
        -----------------
        fit : stanfit4model object
            Sampled model object. Executing "fit" prints out the parameter
            details such as r_hat, diverging samples, as well as cursory
            statistics.
        
        samples: pandas DataFrame
            Tidy pandas DataFrame of sampling output. Each row is an individual
            sample.
        
        stats:  pandas DataFrame
            Tidy pandas DataFrame with summary statistics for each parameter.
        """ 
        if self.stan_model == None:
            self.stan_model == stan_model
        if (self.stan_model == None) & (stan_model == None):
            raise RuntimeError('Stan model not given!')
        # Determine if data has been extracted. If not, do it
        if self.extract == False:
            print('Data not yet defined. Extracting data from matfile...')
            f_looped, dwell, fates = self.extract_data(bead_idx=False)
            print('finished!')
        else:
            f_looped = self.f_looped
            dwell = self.dwell
            fates = self.fates
        
        # Check if the bead IDX is present for cutting fates. If so, reprocess
        if self.bead_idx == True:
            fates = self.cut_beads(bead_idx=False)

        # Define the data dictionary
        data_dict = {'J':int(dwell['replicate'].max()) + 1, 'N':len(dwell),
                    'idx':dwell['replicate'].values, 
                    'total_frames':f_looped['total_frames'],
                    'looped_frames':f_looped['looped_frames'],
                    'n_beads': fates['n_beads'], 'n_cuts':fates['n_cuts'],
                    'dwell_time':dwell['dwell_time_ms']}

        # Load the specified Stan model using the vdj.bayes class
        model = vdj.bayes.StanModel(self.stan_model, data_dict=data_dict, 
                                    force_compile=force_compile)
        fit, samples = model.sample(iter=iter, **sampler_kwargs)
        stats = model.summary()
        return [fit, samples, stats]
        

#%%
mat = ProcessTPM(fname='../../data/analysis_280_analyzed 2.mat',
                 stan_model='../stan/pooled_model.stan')
f_looped, dwell, fates = mat.extract_data()

# %%
fit, samples, stats = mat.run_inference(force_compile=False, iter=5000,
                            sampler_kwargs=dict(control=dict(adapt_delta=0.99)))


                       
#%% [markdown]
# It seems like this function works pretty well. Let's now try to look at the
# inference and see how well it agrees with the data.

#%%
# compute the dwell time ECDFs
fig, ax = plt.subplots(1, 1, figsize=(6,4), dpi=125)
dwell['ecdf'] = dwell.groupby('replicate')['dwell_time_ms'].transform(lambda x: x.rank(method='first') / len(x))

# Plot the replicate distributions
i = 0
for g, d in dwell.groupby('replicate'):
    if i == 0:
        label = 'replicate'
        i += 1
    else:
        label = '__nolegend__'
    d.sort_values('dwell_time_ms', inplace=True)
    ax.step(d['dwell_time_ms'], d['ecdf'], lw=1, alpha=0.2, 
            label=label, color='slategray')

# Plot the pooled distribution
x = np.sort(dwell['dwell_time_ms'].values)
y = np.linspace(0, 1, len(x))
ax.step(x, y, '-', lw=1, color='dodgerblue', label='pooled data')

# Plot the hyperparameter CDF 
r = vdj.stats.compute_hpd(samples['r_cut'] + samples['k_unloop'].values, 0.95)

cdf_low = scipy.stats.expon(loc=0, scale=1/r[0]).cdf(x)
cdf_high = scipy.stats.expon(loc=0, scale=1/r[1]).cdf(x)
ax.fill_between(x, cdf_low, cdf_high, color='tomato', alpha=0.5)
ax.set_xlabel('dwell time [ms]')
ax.set_ylabel('cumulative distribution')
# ax.set_xscale('symlog', linthreshx=1E1)
plt.tight_layout()

#%%
_ = plt.hist(dwell['dwell_time_ms'], bins=500, density=True)
low = r[1] * np.exp(-r[1] * x)
high = r[0] * np.exp(-r[0] * x)
_ = plt.fill_between(x, low, high, alpha=0.4, color='tomato',
zorder=1000)

#%%
import corner
_ = corner.corner(samples[['r_cut', 'k_loop', 'k_unloop']])


#%%
rep = 9 
cdf_low = scipy.stats.expon(loc=0, scale=1/r[0]).cdf(x)
cdf_high = scipy.stats.expon(loc=0, scale=1/r[1]).cdf(x)
plt.fill_between(x, cdf_low, cdf_high, color='tomato')
_rep = dwell[dwell['replicate']==rep]
_rep.sort_values('dwell_time_ms', inplace=True)
plt.step(_rep['dwell_time_ms'], _rep['ecdf'])

#%%
