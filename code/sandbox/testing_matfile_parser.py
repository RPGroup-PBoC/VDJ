# -*- coding: utf-8 -*-
# %%[markdown]
# # Testing Input/Output Methods of TPM Matfiles
#%%
import numpy as np
import pandas as pd
import scipy.io
import vdj.io
import vdj.viz
import vdj.bayes
import matplotlib.pyplot as plt
from functools import reduce
vdj.viz.plotting_style()

# %%[markdown]
# In this notebook, we we test the matfile input-output reading function to make
# sure that everything is beging computed correctly

cd code/sandbox

#%%
_mat = scipy.io.loadmat('../../data/analysis_280_analyzed.mat',
matlab_compatible=True)
_mat.keys()
# %%
# Definitely need to define a class to do all this stuff. 
class ProcessTPM(object):
    def __init__(self, fname=None, matfile=None, framerate=30):
        if (self.fname == None) & (self.matfile != None):
            self.file = matfile
        elif (self.fname != None) & (self.matfile == None):
            self.file = scipy.io.loadmat(matfile)
        else:
            raise RuntimeError("Either `fname` or `matfile` must be specified, but not both!")

        if self.framerate <= 0:
            raise ValueError('framerate must be greater than zero!')
        self.framerate = framerate # in ms/frame

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
        for i, rep in enumerate(mat['ontime_comp'][0]):
            # Find only the times where loopstate was observed
            dwell = rep[rep > 0]
            _df = pd.DataFrame([])
            _df['dwell_time_ms'] = dwell / self.framerate
            _df['replicate'] = i + 1 # indexing replicates by 1 
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
        for i, rep in enumerate(_mat['statetrace_comp'][0]):
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
                            'total_time_ms': total_time * self.framerate,
                            'looped_time_ms': np.sum(rep==2) * self.framerate,
                            'n_beads': n_beads,
                            'replicate': i + 1}, ignore_index=True)
 
        # Make the appropriate entries integers
        df['n_beads'] = df['n_beads'].values.astype(int)
        df['replicate'] = df['replicate'].values.astype(int)
        self.f_looped = df.reset_index() 
        return df
 
    def cut_beads(self, bead_idx=True):
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
            n_beads = len(mat['ontime_comp'][0][i])
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
        else:
            df['n_beads'] = df['n_beads'].values.astype(int)
            df['n_cuts'] = df['n_cuts'].values.astype(int)

        df['replicate'] = df['replicate'].values.astype(int)
        self.fates = df.reset_index()
        return df
 
    def extract_data(self, bead_idx=True):
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
            dwell = self.dwell_time()
        if type(self.fates) == bool:
            fates = self.cut_beads(bead_idx=bead_idx)
        else:
            fates = self.fates
        self.extract = True
        return [f_looped, dwell, fates]
 

    # ##########################################################################
    # INFERENCE 
    # ##########################################################################
    def run_inference(self, model_file=None, run_ppc=False, verbose=True, 
                     **kwargs):

        # Load the specified Stan model using the vdj.bayes class
        model = vdj.bayes.StanModel(model_file, **kwargs)
        
        # Define the data dictionary.
        if self.extract == False:
            f_looped, dwell, fates = self.extract
        data_dict = {}
        
        return None
        

#%%
mat = ParseTPM('../../data/analysis_280_analyzed.mat')

# %%
dwell = mat.dwell_time()
cut = mat.cut_beads(bead_idx=False)
f_looped = mat.fraction_looped()



#%%
out = mat.extract()

p

#%%
pd.merge(out[1], out[2], on='replicate')


#%%
pd.concat(out)

#%%
reduce(lambda left, right: pd.merge(left, right, on=['replicate', 'bead_id']))