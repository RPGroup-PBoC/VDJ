# -*- coding: utf-8 -*_
"""
Read and Compile Data Sets
--------------------------------------------------------------------------------
Author(s): Soichi Hirokawa and Griffin Chure
Last Modified: September 25, 2019
License: MIT

Description
--------------------------------------------------------------------------------
This script reads in the `.mat` files generated from analyzing images of
tethered beads from a TPM experiment. It scrapes the files for the relevant
information, including feature statistics, the date and replicate number, as
well as the mutant and condition information. This script is designed to be
run from the `code/processing` directory and relies on accessing the `.mat`
files from the `data` folder from a relative path.

Execution
--------------------------------------------------------------------------------
This script can be run from the command line as follows:

```
$>  pwd
    vdj_recombination/code/processing
$> python compile_data.py
    PRocessing .mat TPM files... [///////////] 100%
```

"""
import numpy as np
import pandas as pd
from tqdm import tqdm
import vdj.io
import glob


# Glob the file names of all .mat files in data
files = glob.glob('../../data/*.mat')

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
    # V10-95 and V10-96 have the same sequence
    df = df.replace(to_replace='V10-95', value='V10-96')
    df.to_csv(f'../../data/{file}', index=False)

