#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import seaborn as sns
import vdj.io 
import vdj.viz
vdj.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/compiled_looping_events.csv')
data = data[(data['salt']=='Mg') & (data['hmgb1']==80)]

cut_data = pd.read_csv('../../data/pooled_cutting_probability.csv')
cut_data = cut_data[(cut_data['hmgb1'] == 80) & (cut_data['salt']=='Mg')]

cut_posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv')
cut_posts = cut_posts[(cut_posts['hmgb1']==80) & (cut_posts['salt']=='Mg')]

#%%
# Consider only the single mutants. 
points = cut_data[cut_data['n_muts']==1].copy()

# Indicate their x position based on the difference. 
ref =  vdj.io.endogenous_seqs()['WT12rss']
ref_seq = ref[0]
ref_idx = ref[1]
for g in points['mutant'].unique():
    seq = vdj.io.mutation_parser(g)
    loc = np.argmax(ref_idx != seq['seq_idx'])
    points.loc[points['mutant']==g, 'position'] = loc

#%%
wt_12rss
#%%
