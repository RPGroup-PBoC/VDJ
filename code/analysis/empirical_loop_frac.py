# -*- coding:utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bokeh.io
import bokeh.plotting

# Import compiled data
loop_df = pd.read_csv('../../data/compiled_looping_fraction.csv')
#%%
# Create mean looping fraction for each replicate
loop_df['fraction beads looped'] = loop_df['n_loops'] / loop_df['n_beads']
#%%
# Compile empirical mean and standard deviation of looping fraction for each mutant
mutant_names = loop_df['mutant'].unique()
compiled_loop_df = pd.DataFrame()
for mutant in mutant_names:
    mean_fraction = np.mean(loop_df[loop_df['mutant']==mutant]['fraction beads looped'])
    std_fraction = np.std(loop_df[loop_df['mutant']==mutant]['fraction beads looped'])
    compiled_loop_df = compiled_loop_df.append({'mutant' : mutant,
                    'mean fraction looped' : mean_fraction,
                    'std fraction looped' : std_fraction},
                    ignore_index=True)
compiled_loop_df = compiled_loop_df[['mutant', 'mean fraction looped', 'std fraction looped']]
# Save to CSV
compiled_loop_df.to_csv('../../data/compiled_bead_looping_fraction.csv',
                    index=False)
#%%
