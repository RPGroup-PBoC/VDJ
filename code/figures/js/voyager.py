# -*- coding: utf-8 -*-
"""
Builds a Bokeh appelet for exploring the data and visualizing the inference
statistics.
"""
import numpy as np
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, IndexFilter
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.embed import components

# Load the datasets
dwell_times = pd.read_csv('../../../data/compiled_dwell_times.csv')
f_looped = pd.read_csv('../../../data/compiled_looping_fraction.csv')
fates = pd.read_csv('../../../data/compiled_cutting_events.csv')



