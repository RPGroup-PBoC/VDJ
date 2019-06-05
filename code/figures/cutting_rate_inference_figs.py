# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.models
import bokeh.layouts
import scipy.stats
import vdj.io
bokeh.io.output_notebook()

# Load the statistics
data = pd.read_csv('../../data/compiled_dwell_times.csv')
stats = pd.read_csv('../../data/pooled_cutting_rate.csv')


# %%
# Set up the figures
dwell_ax = bokeh.plotting.figure(width=400, height=400, 
                                x_axis_label='dwell time [s]',
                                y_axis_label='cumulative distribution')

rate_ax = bokeh.plotting.figure(width=600, height=200, 
                               x_range=[1, 29], y_range=[1,4],
                               x_axis_label='reference sequence', 
                               y_axis_label='mutation')
rel_rate_ax = bokeh.plotting.figure(width=600, height=200,
                               x_range=[1, 29], y_range=[1, 4],
                               x_axis_label='reference sequence',
                               y_axis_label='mutation')



# Define the layout
col = bokeh.layouts.column(rate_ax, rel_rate_ax)
layout = bokeh.layouts.gridplot([[dwell_ax, col]])
bokeh.io.show(layout)


#%%


#%%
