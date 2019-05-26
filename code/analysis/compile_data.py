# -*- coding: utf-8 -*_
"""
Reads in a set of compiled data for all mutants, runs processing, and saves
relevant information to .csv. This does *not* perform inference on the files --
only generates the csv files. 
"""
import numpy as np
import pandas as pd
from tqdm import tqdm
import vdj.io
import glob

# Glob the file names of all .mat files in data
files = glob.glob('../../data/*.mat')

# Instantiate the empty lists to store dataframes
dwell_dfs, f_looped_dfs, fates_dfs = [], [], []
for f in tqdm(files, desc='Processing .mat TPM files...'):
    f_looped, dwell, fates = vdj.io.ProcessTPM(f).extract_data()
    dwell_dfs.append(dwell)
    f_looped_dfs.append(f_looped)
    fates_dfs.append(fates)

# Concatenate the data and save to CSVs
dwell = pd.concat(dwell_dfs)
dwell.to_csv('../../data/compiled_dwell_times.csv', index=False)
f_looped = pd.concat(f_looped_dfs)
f_looped.to_csv('../../data/compiled_looping_fraction.csv', index=False)
fates = pd.concat(fates_dfs)
fates.to_csv('../../data/compiled_cutting_events.csv', index=False)
