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
files = glob.glob('../../data/*.mat')

# Instantiate the empty lists to store dataframes
dwell_dfs, f_looped_dfs, fates_dfs = [], [], []
samples_dfs, stats_dfs = [], []
for f in tqdm(files, desc='Processing .mat TPM files...'):
    # Load and extract data
    tpm = vdj.io.ProcessTPM(f)
    f_looped, dwell, fates = tpm.extract_data()

    # Perform the hierarchical inference.
    _, samples, stats = tpm.run_inference(
        stan_model='../stan/hierarchical_model.stan',
        iter=5000, sampler_kwargs=dict(n_jobs=1, control=dict(adapt_delta=0.9)))

    # Append the generated dataframes to the data lists
    dwell_dfs.append(dwell)
    f_looped_dfs.append(f_looped)
    fates_dfs.append(fates)

    # Append the inference data
    samples_dfs.append(samples)
    stats_dfs.append(stats)

# Concatenate the data and save to CSVs
dwell = pd.concat(dwell_dfs)
dwell.to_csv('../../data/compiled_dwell_times.csv', index=False)
f_looped = pd.concat(f_looped_dfs)
f_looped.to_csv('../../data/compiled_looping_fraction.csv', index=False)
fates = pd.concat(fates_dfs)
fates.to_csv('../../data/compiled_cutting_events.csv', index=False)
samples = pd.concat(samples_dfs)
samples.to_csv('../../data/compiled_inference_samples.csv', index=False)
stats = pd.concat(stats_dfs)
stats.to_csv('../../data/compiled_inference_statistics.csv', index=False)
