# -*- coding: utf-8 -*_
"""
Reads in a set of compiled data for all mutants, runs processing, and saves
relevant information to .csv. 
"""
import numpy as np
import pandas as pd
from tqdm import tqdm
import vdj.io
import glob
import imp
imp.reload(vdj.io)

# Glob the file names of all .mat files in data
files = glob.glob('../../data/mat_files/*.mat')

# Instantiate the empty lists to store dataframes
dwell_dfs, f_looped_dfs, fates_dfs, events_dfs = [], [], [], []
dfs = [f_looped_dfs, dwell_dfs, fates_dfs, events_dfs]
for f in tqdm(files, desc='Processing .mat TPM files...'):
    # Load and extract data
    tpm = vdj.io.ProcessTPM(f)
    out = tpm.extract_data()

    # Add divalent salt information
    if '_Ca_' not in f.split('/')[-1]:
        salt = 'Mg'
    else:
        salt = 'Ca'
    for d in out:
        d['salt'] = salt

    # Add HMGB1 concentration information
    if '_HMGB1_' not in f.split('/')[-1]:
        hmgb1 = 80
    else:
        hmgb1 = int(f.split('_')[2][:-2])
    for d in out:
        d['hmgb1'] = hmgb1

    # Append the generated dataframes to the data lists
    for i, d in enumerate(dfs):
        d.append(out[i])

# Concatenate the data and save to CSVs
fnames = ['compiled_looping_fraction.csv', 'compiled_dwell_times.csv',
         'compiled_bead_fates.csv', 'compiled_looping_events.csv']
for file, d in zip(fnames, dfs):
    df = pd.concat(d)
    # V4-55 RSS and 12SpacC1A have the same sequence
    df = df.replace(to_replace='V4-55', value='12SpacC1A')
    # V10-95 and V10-96 have the same sequence
    df = df.replace(to_replace='V10-95', value='V10-96')
    df.to_csv(f'../../data/{file}', index=False)


