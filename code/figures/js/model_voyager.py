# -*- coding: utf-8 -*-
"""
Builds a Bokeh appelet for exploring the model and 
statistics.
"""
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, GroupFilter 
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider
from bokeh.models.glyphs import HBar
from bokeh.embed import components
bokeh.plotting.output_file('./model_voyager.html')
import vdj.io
bokeh.io.output_notebook()


# ##############################################################################
# INSTANTIATE STATIC ELEMENTS
# ##############################################################################
# Define the sliders. 
rcut_slider = Slider(title='log10 (cutting rate [s^-1])', start=-10, end=0, 
                    step=0.1, value=-3)
kloop_slider = Slider(title='log10 (looping rate [s^-1])', start=-10, end=0, 
                    step=0.1, value=-5)
kunloop_slider = Slider(title='log10 (unlooping rate [s^-1])', start=-10, end=0, 
                    step=0.1, value=-2)

# Set up the plotting canvas
pdf_ax = bokeh.plotting.figure(width=300, height=300, 
                                x_axis_label='paired complex dwell time [s]',
                                y_axis_label='probability density')
cdf_ax = bokeh.plotting.figure(width=300, height=300,
                                  x_axis_label='paired complex dwell time [s]',
                                  y_axis_label='cumulitative density')
prob_ax = bokeh.plotting.figure(width=300, height=200, y_axis_label='probability')

# #############################################################################
# DEFINE THE INITIAL CONDITIONS
# #############################################################################
# Dwell time distributions
time_range = np.logspace(0, 4, 500)
r = rcut_slider.value + kunloop_slider.value
dwell_pdf = r * np.exp(-r * time_range)
dwell_cdf = 1 - np.exp(-r * time_range)

#Probabilities
rcut = 10**rcut_slider.value
kloop = 10**kloop_slider.value
kunloop = 10**kunloop_slider.value
pcut = [rcut / (rcut + kunloop)]
floop = [kloop / (rcut + kunloop + kloop)]
ploop = [kloop / (kloop + kunloop)]

# ##############################################################################
# DEFINE THE DATA SOURCES
# ##############################################################################
dwell_source = ColumnDataSource({'time':time_range, 'pdf':dwell_pdf, 
                                'cdf':dwell_cdf})
prob_source = ColumnDataSource({'P_loop':ploop, 'P_cut':pcut, 'F_loop':floop})

# ##############################################################################
# DEFINE JS CALLBACK
# ##############################################################################
cb = CustomJS(args={'prob_source':prob_source, 'dwell_source':dwell_source,
                    'rcut_slider':rcut_slider, 'kloop_slider':kloop_slider,
                    'kunloop_slider':kunloop_slider}, 
code = """
// Variable definition
var prob_source = prob_source.data;
var dwell_time = dwell_source.data['time'];
var dwell_cdf = dwell_source.data['cdf'];
var dwell_pdf = dwell_source.data['pdf'];
var r_cut = Math.pow(10, rcut_slider.value);
var k_loop = Math.pow(10, kloop_slider.value);
var k_unloop = Math.pow(10, kunloop_slider.value);
var r = r_cut + k_unloop;


// Compute the quantities of interest
prob_source['p_loop'] = k_loop / (k_loop + k_unloop);
prob_source['p_cut'] = r_cut / (r_cut + k_unloop);
prob_source['f_loop'] = k_loop / (r + k_loop); 

// Compute the dwell time cumulative distribution
for (var i = 0; i < dwell_cdf.length; i++) {
    dwell_pdf[i] = r * Math.exp(-r * dwell_time[i]);
    dwell_cdf[i] = 1 - Math.exp(-r * dwell_time[i]);
} 
dwell_pdf = sum(dwell_pdf);
prob_source.change.emit();
dwell_source.change.emit(); 
""")

# Assign callbacks to sliders. 
for s in [rcut_slider, kloop_slider, kunloop_slider]:
    s.callback = cb
    s.js_on_change('value', cb)

# ##############################################################################
# POPULATE AXES
# ##############################################################################
pdf_ax.line(x='time', y='pdf', line_width=2, color='tomato', 
               source=dwell_source)
cdf_ax.step(x='time', y='cdf', line_width=2, color='tomato', 
               source=dwell_source)



box = widgetbox([rcut_slider, kloop_slider, kunloop_slider])
layout = bokeh.layouts.gridplot([[box, prob_ax], [pdf_ax, cdf_ax]])
# #############################################################################
#  THEME DETAILS
# #############################################################################
theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#f5e3b3',
                'outline_line_color': '#000000',
            },
            'Axis': {
            'axis_line_color': "black",
            'major_tick_out': 7,
            'major_tick_line_width': 0.75,
            'major_tick_line_color': "black",
            'minor_tick_line_color': "black",
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': None,
            },
            'Legend': {
                'background_fill_color': '#f5e3b3',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica',
                'offset': 2,
            }}}

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(layout)



