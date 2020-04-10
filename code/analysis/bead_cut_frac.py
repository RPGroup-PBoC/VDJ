# -*- coding: utf-8 -*-
#%% 
import numpy as np
import pandas as pd
from scipy.special import gammaln, logsumexp
import vdj.io


# Load the fates of beads and number of beads
cut_data = pd.read_csv('../../data/compiled_bead_fates.csv')
bead_data = pd.read_csv('../../data/compiled_looping_fraction.csv')

# Get the mutant info
mut_info = {m:vdj.io.mutation_parser(m) for m in cut_data['mutant'].unique()}

# Compute the pooled statistics
pooled_cut = cut_data.groupby(['mutant', 'salt', 'hmgb1']).agg(('sum')).reset_index()
pooled_cut = pooled_cut[['mutant', 'salt', 'hmgb1', 'n_beads', 'n_cuts']]
pooled_cut = pooled_cut.rename(columns={'n_beads' : 'n_loops'})

pooled_bead = bead_data.groupby(['mutant', 'salt', 'hmgb1']).agg(('sum')).reset_index()
pooled_bead = pooled_bead[['mutant', 'salt', 'hmgb1', 'n_beads']]

pooled = pd.merge(pooled_bead, pooled_cut, on=['mutant','salt','hmgb1'])

pooled['mode'] = pooled['n_cuts'].values / pooled['n_beads']
pooled['std'] = np.sqrt((pooled['n_cuts'].values * (pooled['n_beads'] -\
                 pooled['n_cuts'])) / pooled['n_beads'].values**3)
for m, seq in mut_info.items():
    pooled.loc[pooled['mutant']==m, 'n_muts'] = seq['n_muts']
pooled.to_csv('../../data/pooled_bead_cutting_probability.csv', index=False)

# Compute the hierarchical stats
hier = pooled.copy()
hier['mode'] = hier['n_cuts'].values / hier['n_beads']
hier['std'] = (hier['n_cuts'].values * (hier['n_beads'] -\
             hier['n_cuts'].values)) / hier['n_beads']**3
for m, seq in mut_info.items():
    hier.loc[hier['mutant']==m, 'n_muts'] = seq['n_muts']
hier.to_csv('../../data/independent_bead_cutting_probability.csv', index=False)


