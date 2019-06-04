import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.figure
import vdj.io
import vdj.stats
bokeh.io.output_notebook()

#%% 
# Load the sampling information
cut_stats = pd.read_csv('../../data/cutting_probability_summary.csv')

# Get the reference sequence
seqs = vdj.io.endogenous_seqs()
ref_seq  = seqs['reference'][0]
ref_int = seqs['reference'][1]

# isolate the point mutants
point_muts = cut_stats[cut_stats['n_muts'] == 1]

# Define the nucleotide position
for m in point_muts['mutant'].unique():
    # Parse the mutation 
    mut = vdj.io.mutation_parser(m)

    # Find the index where they don't agree
    ind = np.argmax()


# %%
ax = bokeh.plotting.figure(width=900, height=300, x_axis_label='probability',
                          y_axis_label='reference sequence')




#%%



