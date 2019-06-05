"""
io.py

This submodule contains an array of utilities for simple file input/output and
other file parsing utilities 
"""
import numpy as np
import pandas as pd
from .bayes import StanModel
import scipy.io
import os

def nucleotide_idx():
    """
    Returns the dictionary linking base identity to an integer
    """
    return {'A':0, 'C':1, 'G':2, 'T':3, 0:'A', 1:'C', 2:'G', 3:'T'}

def endogenous_seqs():
    """
    Returns a dictionary of the sequence identity for the reference 12RSS
    sequence and useful mutations. 
    """
    conv = nucleotide_idx()
    _seqs = {'reference': 'CACAGTGCTACAGACTGGAACAAAAACC',
            'WT12rss': 'CACAGTGCTACAGACTGGAACAAAAACC',
            'V1-135':  'CACAGTGATTCAGACCCGAACAAAAACT',
            'V9-120': 'CACAGTGATACAAATCATAACATAAACC',
            'V8-18': 'CACAGAGCTTCAGCTGCCTACACAAACC',
            'V6-17': 'CACAGTGCTTCAGCCTCCTACACAAACC',
            'V10-96': 'CACAATGATATAAGTCATAACATAAACC',
            'V10-95': 'CACAATGATATAAGTCATAACATAAACC',
            'V19-93': 'CACAGTGATACAAATCATAACAAAAACC',
            'V4-55': 'CACAGTGATACAGACTGGAACAAAAACC',
            'V5-43': 'CACAGTGATGCAGACCATAGCAAAAATC',
            'V6-15': 'CACAGTACTTCAGCCTCCTACATAAACC',
            'DFL161': 'CACAGTGCTATATCCATCAGCAAAAACC',
            'DFL1613': 'CACAGTAGTAGATCCCTTCACAAAAAGC'}

    # Add the integer conversions
    ref_seq = _seqs['reference']
    seqs = {m: [seq, np.array([conv[a] for a in seq]),  
        np.sum(np.array(list(ref_seq)) != np.array(list(seq)))] for m, seq in _seqs.items()}

    return seqs

def mutation_parser(mut_id):
    """
    Takes a mutation identifier and returns the sequence in both base and
    integer representation

    Parameters
    ----------
    mut_id: str
        String of the mutant identifier with the region ('Hept', 'Spac', or
        'Non') and mutation(s) such as A4T.
    
    Returns
    -------
    mut_dict: dict
        Dictionary with the supplied mut_id, new 12RSS sequence, and integer
        representation of the new sequence.
    """
    # Load the raw sequences and nt conversion id's
    conversion = nucleotide_idx()
    seqs = endogenous_seqs()
    ref = seqs['reference'][0]

    # Determine the region which is mutated
    if 'hept' in mut_id.lower():
        seq = ref[:7]
        region = 'hept'
    elif 'spac' in mut_id.lower():
        seq = ref[7:19]
        region =  'spac'
    elif 'non' in mut_id.lower():
        seq = ref[19:]
        region = 'non'
    else:
        try:
            new_seq = seqs[mut_id]
            return {'seq':new_seq[0], 'seq_idx':new_seq[1], 'n_muts':new_seq[2]}
        except:
            raise ValueError("Mutation ID not in proper format of 12(region)(old)(pos)(new) or name not recognized.")
            
    # Force the string to uppercae
    loc = mut_id.lower().split(region)[-1].upper()

    # Get the indices of base positions
    base_id = np.array([i for i, b in enumerate(list(loc)) if b in 'ATCG']).astype(int)

    # Split the indices into pairs
    if (len(base_id) / 2) == 1:
        muts = [loc]
    else:
        inds = [(start, end) for start, end in zip(base_id[::2], base_id[1::2])]
        muts = [loc[ind[0]:ind[1] + 1] for ind in inds]

    new_region = list(seq)
    for m in muts:
        pos = int(m[1:-1])
        if seq[pos - 1] != m[0].upper():
            raise ValueError(f'Base at position {m[1]} is not {m[0].upper()}! Double check the sequence.')
        else:
            new_region[pos - 1] = m[-1]

    # Replace the region and compute the integer representation
    new_seq = ref.replace(seq, ''.join(new_region).upper())
    new_seq_idx = np.array([conversion[a] for a in new_seq]) 
    return {'seq':new_seq, 'seq_idx':new_seq_idx, 'n_muts':len(muts)}


class ProcessTPM(object):
    def __init__(self, fname=None, framerate=30, 
                 stan_model=None):
        if framerate <= 0:
            raise ValueError('framerate must be greater than zero!')
        if type(fname) != str:
            raise TypeError('File must be a string')

        # Load the matfile
        self.file = scipy.io.loadmat(fname)

        # Define variables
        self.fps = framerate # n ms/frame
        self.stan_model = stan_model
        self.mut_id = fname.split('/')[-1].split('_')[0]
        self.seq = mutation_parser(self.mut_id)['seq']

       # Generate a dictionary of dates and add to self.
        try:
            experiment_names = self.file['alllacnames'][0]
            dates = {}
            for i, nom in enumerate(experiment_names):
                dates[i+1] = nom[0].split('_')[0]
            self.dates = dates 
        except KeyError:
            self.dates = {i+1:'Date Unknown' for i in range(1000)}

        # Add extra "dates unknown"
        for i in range(1000):
            self.dates[i + np.max(list(self.dates.keys()))] = 'Date Unknown'

        # Initialize state class state variables
        self.drop_replicate = None # If a replicate is missing
        self.extract = False # If data has been extracted yet from matfile

        # Initialize the extracted quantities
        self.dwell = False
        self.f_looped = False
        self.fates = False

    # ##########################################################################
    # DATA EXTRACTION
    # ##########################################################################
    def dwell_time(self):
        """
        Extracts the dwell times for each replicate.
        """
        mat = self.file
        dfs = [] # empty list to append coming data frames

        for i, rep in enumerate(mat['loops'][0]):
            # Find only the times where loopstate was observed
            dwell = rep[rep > 0]

            # Isolate the cut dwell times for the replicate 
            cuts = mat['pccut'][0][i]
            cuts = cuts[cuts > 0]
            
            # If there were no looped states, increment the counter and move on
            _dfs = []
            if len(dwell) == 0:
                dwell = np.array([0])
                cut = np.array([0])
            else:
                # Find out which dwell times lead to a cutting event
                if len(cuts) == 0:
                    _df = pd.DataFrame([])
                    _df['dwell_time_s'] = dwell / self.fps
                    _df['cut'] = 0
                    _df['replicate'] = i + 1
                    _dfs.append(_df)
                else:
                    cut_dwells = dwell[cuts==dwell]
                    unloop_dwells = dwell[cuts!=dwell]
                    self.cut_dwells = cut_dwells
                    self.unloop_dwells = unloop_dwells


                    for c in cut_dwells:
                        _df = pd.DataFrame([])
                        _df['dwell_time_s'] = c / self.fps, 
                        _df['cut']  = 1
                        _df['replicate'] = i + 1
                        _dfs.append(_df)
                    if len(unloop_dwells) > 0:
                        for u in unloop_dwells:
                            _df = pd.DataFrame([])
                            _df['dwell_time_s'] = u / self.fps
                            _df['cut'] = 0
                            _df['replicate'] = i + 1
                            _dfs.append(_df)
                dfs.append(pd.concat(_dfs, sort=False))
        df = pd.concat(dfs, sort=False).reset_index()
 
        # Make the appropriate entries integers
        df['mutant'] = self.mut_id
        df['replicate'] = df['replicate'].values
        df['date'] = self.dates[i]
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
            _dwells = mat['loops'][0][i]
            n_loops = len(_dwells[_dwells > 0])
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
                            'n_loops': n_loops,
                            'replicate': i + 1}, ignore_index=True)
 
        # Make the appropriate entries integers
        df['mutant'] = self.mut_id
        df['n_beads'] = df['n_beads'].values.astype(int)
        df['replicate'] = df['replicate'].values.astype(int)
        df['total_frames'] = df['total_frames'].values.astype(int)
        df['looped_frames'] = df['looped_frames'].values.astype(int)
        df['n_loops'] = df['n_loops'].values.astype(int)

        df['date'] = self.dates[i]
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
                                'replicate': i+1}, ignore_index=True)
                
        # Make the appropriate entries integers
        if bead_idx == True:
            df['bead_id'] = df['bead_id'].values.astype(int)
            self.bead_idx = True
        else:
            df['n_beads'] = df['n_beads'].values.astype(int)
            df['n_cuts'] = df['n_cuts'].values.astype(int)
            self.bead_idx = False

        df['replicate'] = df['replicate'].values.astype(int)
        df['mutant'] = self.mut_id
        df['date'] = self.dates[i]
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
 

