"""
Generates SI Figure for bulk DNA tether loss with data on looping fraction
and cutting probability.
"""
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib.lines as lines
import vdj.viz
import vdj.io
vdj.viz.plotting_style()

# Load the data with long-form looping events and restrict to relevant sets.
data = pd.read_csv('../../data/compiled_loop_freq_bs.csv',       
                   comment='#')
counts = data[(data['salt']=='Mg') & (data['hmgb1']==80) & 
                (data['percentile']==95.0) & (data['mutant']!='12CodC6A')]

# Load all cutting probability estimates taking Gaussian approximation.
cut_data = pd.read_csv('../../data/pooled_cutting_probability.csv', comment='#')
cut_data = cut_data[(cut_data['hmgb1'] == 80) & (cut_data['salt']=='Mg') & 
                (cut_data['mutant']!='12CodC6A')]

# Load all bead cleavage probability estimates taking Gaussian approximation.
bead_data = pd.read_csv('../../data/pooled_bead_cutting_probability.csv', comment='#')
bead_data = bead_data[(bead_data['hmgb1'] == 80) & (bead_data['salt']=='Mg') & 
                        (bead_data['mutant']!='12CodC6A')]


