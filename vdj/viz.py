"""
viz.py

Module for calling standard plotting functions and setting the styles.
"""
import matplotlib.pyplot as plt
import seaborn as sns
from . import io
from bokeh.io import curdoc
from bokeh.themes import Theme

def plotting_style():
    """
    Sets matplotlibs plotting style to the publication default. It returns a
    list of the preferred colors.
    """
    rc = {'axes.facecolor': '#f5e3b3',
          'axes.grid': False,
          'axes.frameon': True,
          'ytick.direction': 'out',
          'xtick.direction': 'out',
          'xtick.gridOn': True,
          'ytick.gridOn': True,
          'ytick.major.width':5,
          'xtick.major.width':5,
          'ytick.major.size': 5,
          'xtick.major.size': 5,
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.facecolor': '#FFEDCE',
          'figure.dpi': 150,
          'xtick.color': 'k',
          'ytick.color': 'k'}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('dark', rc=rc)

def plotting_style_bokeh():
    """
    Sets bokeh plotting style to the publication default. It returns a list
    of the preferred colors.
    """
    curdoc().theme = Theme(json={'attrs' : {
        
        # apply defaults to Figure properties
        'Figure': {
            'toolbar_location': None,
            'outline_line_color': None,
            'min_border_right': 10,
        },

        # apply defaults to Axis properties
        'Axis': {
            'major_tick_in': None,
            'minor_tick_in': None,
        }
    }})


def generate_matrix(mutations, data_values):
    """
    Generates a grid of values the size of the 12RSS sequence

    Parameters
    ----------
    mutation_array : 
    """
    seqs = io.endogenous_seqs()
    consensus = seqs['consensus']
    return None
