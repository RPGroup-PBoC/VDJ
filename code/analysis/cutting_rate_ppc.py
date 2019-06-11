# -*- coding: utf-8 -*-
# %% [markdown]
# Prior Predictive Checks for Cutting Rate Inference
# %%
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import scipy.stats as st
bokeh.io.output_notebook()


# %% [markdown]
# I've noticed in playing around with inference of the cutting rate that the
# choice of prior significantly influences resulting postieror. It seems like
# the trend -- that most of the mutations *increase* the cutting rate holds, the
# exact values vary. In this notebook, we explore the prior choices for $\tau$
# through prior predictive checks to narrow down on the best choice for the
# prior.
#
# As a reminder, our statistical model for inference of the cutting rate $\tau$
# comes from Bayes' theorem as
# $$ g(\tau\,\vert\, \{t\}, \{T\}) \propto f(\{t\}\,\vert\, \tau)f(\{T\}, \vert\, \tau)g(\tau)$$
# the likelihoods can be written as 
# $$f(\{t\}\,\vert\, \tau) = {1 \over \tau^N}e^{-N\bar{t}\over\tau}$$
# and 
# $$f(\{T\}\,\vert\,\tau) = e^{-M\bar{T}\over \tau}$$
# Where $N$ and $M$ correspond to the number of cuts and unlooping events,
# respectively. 
#
# The choice of prior is the hard part. Let's explore a few of them. 
#%% [markdown]

# ## Defining the log posterior.

# As we don't have to actually sample this and can just plot the analytical
# solution over the range of $\tau$ we care about, we'll start by defining a
# function to compute the log posterior.

# %%
def log_like(t, T, tau):
    N = len(t)
    M = len(T)
    return -N * np.log(tau) * (-N * np.mean(t) -M * np.mean(T)) / tau

# %% [markdown]

# We can also define some constants, such as the range of tau, over which we 
# %% [markdown]
# ## Choice 1: Half-Normal
