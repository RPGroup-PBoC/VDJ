#%%
import numpy as np
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.events import Tap
from bokeh.models import * 
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select 
from bokeh.embed import components
import vdj.io
import vdj.stats
bokeh.io.output_notebook()
bokeh.plotting.output_file('./pcut_model_explorer.html')
# %%
# Define the sliders
n_loops_slider = Slider(start=1, end=500, step=1, value=10,
                title='number of looping events')
n_cuts_slider = Slider(start=1, end=n_loops_slider.value, step=1, value=1,
                title='number of cutting events')

# Define the parameters to evaluate
prob_range = np.linspace(0, 1, 500)
post_source = ColumnDataSource({'p':prob_range, 'posterior': np.zeros_like(prob_range)})

# Define the callback args. 
args = {'post_source':post_source, 
        'n_loops_slider':n_loops_slider, 
        'n_cuts_slider':n_cuts_slider}

# Define the callback
cb = CustomJS(args=args, code="""
// Define a gamma function using the Lanczos approximation
// See https://github.com/substack/gamma.js for node module
 var g_ln = 607/128;
 var p_ln = [0.99999999999999709182,57.156235665862923517,-59.597960355475491248,
            14.136097974741747174,-0.49191381609762019978,0.33994649984811888699e-4,
            0.46523628927048575665e-4,-0.98374475304879564677e-4,0.15808870322491248884e-3,
            -0.21026444172410488319e-3,0.21743961811521264320e-3,-0.16431810653676389022e-3,
            0.84418223983852743293e-4, -0.26190838401581408670e-4,0.36899182659531622704e-5];
 // Spouge approximation (suitable for large arguments)
function lngamma(z) {  
      if (z < 0) {return Number('0/0')};
      var x = p_ln[0];
      for(var i = p_ln.length - 1; i > 0; --i) { 
            x += p_ln[i] / (z + i)
        };
      var t = z + g_ln + 0.5;
      return .5*Math.log(2*Math.PI)+(z+.5)*Math.log(t)-t+Math.log(x)-Math.log(z);
  }
// Define the log likelihood
function logLike(n, N, p) {
    var binomCoeff = lngamma(N + 2) - lngamma(n + 1) - lngamma(N - n + 1)
    var prob = n * Math.log(p) +  (N - n) * Math.log(1 - p)
    return binomCoeff + prob;
}

function logSumExp(vals) {
    var maxVal = Math.max.apply(null, vals);
    var sum_value = 0
    for (var i = 0; i < vals.length; i++) {
        sum_value += Math.exp(vals[i] - maxVal);
    } 
   return maxVal + Math.log(sum_value);
}

// Get the input values
var N = n_loops_slider.value;
var n = n_cuts_slider.value;

// Set the maximum of the cutting slider to the maximum of the loop slider
n_cuts_slider.end = n_loops_slider.value - 1

// Evaluate the posterior
var log_posterior = [];

for (var i = 0; i < post_source.data['p'].length; i++) {
    var prob = post_source.data['p'][i]
    var log_likelihood = logLike(n, N, prob);
    log_posterior[i] = log_likelihood;
}

// Use the logsumexp trick to normalize the posterior
var summed = logSumExp(log_posterior);
var posterior= []
for (var i = 0 ; i < log_posterior.length; i++) { 
   posterior[i] = Math.exp(log_posterior[i] - summed);
}

// Rescale everything to be on the same meaningless scale ([0, 1]).
var maxPosterior = Math.max(...posterior)
var posteriorNorm = []
for (var i = 0; i < post_source.data['p'].length; i++) {
    posteriorNorm[i] = posterior[i] / maxPosterior;
}

// Update the source data.
post_source.data['posterior'] = posteriorNorm;
post_source.change.emit();
"""
)

# Assign the callbacks to the sliders on change. 
n_loops_slider.js_on_change('value', cb)
n_cuts_slider.js_on_change('value', cb)

# Set up the plotting axis.
p = bokeh.plotting.figure(width=600, height=400, x_axis_label='cutting probability',
                        y_axis_label = 'posterior probability', y_range=[0,1],
                        tools=[])
p.line('p', 'posterior', color='dodgerblue', source=post_source)
p.varea('p', 0, 'posterior', source=post_source, fill_color='dodgerblue',
            alpha=0.25)
p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
p.yaxis.minor_tick_line_color = None 
p.yaxis.major_label_text_font_size = '0pt'  # turn off y-axis tick labels

lay = bokeh.layouts.column(n_loops_slider, n_cuts_slider, p)


# Add the theme
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
                'text_font_style': 'bold',
                'align': 'center',
                'text_font': 'Helvetica',

                'offset': 2,
            }}}

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(lay)
# Link the cuts slider to the loops slider to 
# %%
