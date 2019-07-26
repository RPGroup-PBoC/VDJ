# %%[markdown]
# # Analyzing RMS Bead Traces using Bayesian Changepoint Analysis
#%%
import numpy as np
import pandas as pd
import scipy.io
import bokeh.io
import bokeh.plotting
import bokeh.models
import scipy.special
import vdj.viz 
import vdj.bayes
import vdj.stats
import scipy.io
vdj.viz.plotting_style()
bokeh.io.output_notebook()

# %%[markdown]
# In this notebook, we explore the application of Bayesian change point analysis
# to better identify when there are changes in the statistics of RMS
# trajectories. If we can get this to work, we can apply it to teh data in which
# only the 12RSS is present to measure binding kinetics. 

#%%
# step 1 load the data. 
mat = scipy.io.loadmat('../../data/lacdata.mat')
# %%
# Get this into a tidy dataframe -- this should definitely be turned into a
# tested function
FRAME_RATE = 30 # in ms
dfs = []
bead_id = 0
for rep in range(np.shape(mat['lacANAL_RMS'][0][0])[0] - 1):
        samp = mat['lacANAL_RMS'][0][0][rep]
        bead_id += 1
        _df = pd.DataFrame([])
        _df['rms_nm'] = samp
        _df['time_ms'] = np.arange(0, len(samp), 1)
        _df['bead_id'] = bead_id
        dfs.append(_df)
df = pd.concat(dfs)

# Save it to a csv.
# df.to_csv('../../data/consensus_trace_tidy.csv', index=False)

# %% 
# Plot some of the traces
# GOOD_IDX = [99, 100]
p = bokeh.plotting.figure(width=500, height=300, x_axis_label='time [ms]',
                        y_axis_label='RMS')

palette = bokeh.palettes.Category20_20
iter = 0
for g, d in df.groupby(['bead_id']):
    p.line(d['time_ms'], d['rms_nm'], color=palette[iter])
    iter += 1
bokeh.io.show(p)

#%%[markdown]
# Wow that is really, really subtle. Let's try to cut out the trace which has
# rMS of 0
#%%
df = df[df['rms_nm'] > 0]
p = bokeh.plotting.figure(width=500, height=300, x_axis_label='time [ms]',
                        y_axis_label='RMS')

palette = bokeh.palettes.Category20_20
iter = 0
for g, d in df.groupby(['bead_id']):
    if g == 2:
        p.line(d['time_ms'], d['rms_nm'], color=palette[iter])
        iter += 1
bokeh.io.show(p)

#%% 
model = vdj.bayes.StanModel('../stan/changepoint.stan') 

#%%
# Choose one example.
bead = 2
_df = df[df['bead_id']==bead]

data_dict = {'N':len(_df['rms_nm'].values[::100]), 'rms':_df['rms_nm'].values[::100]}
sample = model.sample(data_dict)

#%%

