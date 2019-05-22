# %%[markdown]
# # Analyzing RMS Bead Traces using Bayesian Changepoint Analysis

#%%
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt 
import vdj.viz 
import vdj.bayes
import vdj.stats
import scipy.io
vdj.viz.plotting_style()

# %%
# step 1 load the data. 
mat = scipy.io.loadmat('../../data/lacdata.mat')
# %%
# Get this into a tidy dataframe -- this should definitely be turned into a
# tested function
FRAME_RATE = 30 # in ms
dfs = []
bead_id = 0
for rep in range(np.shape(mat['lacANAL_RMS'])[1] - 1):
    n_beads = np.shape(mat['lacANAL_RMS'][0][rep])[0]
    for i, trace in enumerate(mat['lacANAL_RMS'][0][rep]):
        bead_id += 1
        _df = pd.DataFrame([])
        _df['rms_nm'] = trace
        _df['time_ms'] = np.arange(0, len(trace), 1) * FRAME_RATE
        _df['bead_id'] = bead_id
        _df['replicate'] = rep + 1
        dfs.append(_df)
df = pd.concat(dfs)

# Save it to a csv.
df.to_csv('../../data/consensus_trace_tidy.csv', index=False)



# %% 
# Plot some of the traces
GOOD_IDX = [99, 100]
fig, ax = plt.subplots(1, 1)
sub = df[(df['bead_id'] == 99)]
for g, d in sub.groupby('bead_id'):
    ax.plot(d['time_ms'], d['rms_nm'], '-')


#%% 
model = vdj.bayes.StanModel('../stan/changepoint.stan')

#%%
data_dict = {'N':len(sub), 'rms':sub['rms_nm']}
sample = model.sample(data_dict)

#%%
