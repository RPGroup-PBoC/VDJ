# -*- coding:utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bokeh.io
import bokeh.plotting

# Import compiled data
cut_df = pd.read_csv('../../data/compiled_cutting_events.csv')
#%%
# Create column for mean fraction of loops cut for each replicate
cut_df['cutting_fraction'] = cut_df['n_cuts'] / cut_df['n_beads']
#%%
# Compile empirical mean and standard deviation of cutting fraction for each mutant
mutant_names = cut_df['mutant'].unique()
compiled_cut_df = pd.DataFrame()
for mutant in mutant_names:
    mean_fraction = np.mean(cut_df[cut_df['mutant']==mutant]['cutting_fraction'])
    std_fraction = np.std(cut_df[cut_df['mutant']==mutant]['cutting_fraction'])
    global_fraction = np.sum(cut_df[cut_df['mutant']==mutant]['n_cuts']) / np.sum(cut_df[cut_df['mutant']==mutant]['n_beads'])
    compiled_cut_df = compiled_cut_df.append({'mutant' : mutant,
                    'mean_cuts_per_bead' : mean_fraction,
                    'std_cuts_per_bead' : std_fraction,
                    'global_cuts_per_bead' : global_fraction},
                    ignore_index=True)
compiled_cut_df = compiled_cut_df[['mutant', 'global_cuts_per_bead', 'mean_cuts_per_bead', 'std_cuts_per_bead']]
# Save to CSV
compiled_cut_df.to_csv('../../data/compiled_bead_cutting_fraction.csv',
                    index=False)
#%%
cut_df.to_csv('../../data/compiled_cutting_events.csv')

#%%
