#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import vdj.io
import vdj.viz
vdj.viz.plotting_style()

# Upload data on coding flank-relevant RSSs
cf_muts = ['WT12rss', '12CodC6A', '12SpacC1A', 'V4-55']

loop = pd.read_csv('../../data/compiled_loop_freq_bs.csv')
loop = loop[(loop['mutant'].isin(cf_muts)) & (loop['hmgb1']==80) & (loop['salt']=='Mg')]

dwell = pd.read_csv('../../data/compiled_dwell_times.csv')
dwell = dwell[(dwell['mutant'].isin(cf_muts)) & (dwell['hmgb1']==80) & (dwell['salt']=='Mg')]

cut_posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv')
cut_posts = cut_posts[(cut_posts['mutant'].isin(cf_muts)) & (cut_posts['hmgb1']==80) & (cut_posts['salt']=='Mg')]

#%%


#%%
