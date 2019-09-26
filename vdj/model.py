"""
model.py

Utilities for easy calling of model functions and parameters.
"""
import numpy as np
import scipy.stats as st
import pandas as pd



def unlooping_rate_log_likelihood(tau, t_cut, t_unloop):
    if len(t_cut) == 0:
        n_cuts = 0
        t_cut = np.array([0])
    else:
        n_cuts = len(t_cut)
    if len(t_unloop) == 0:
        n_unloop = 0
        t_unloop = np.array([0])
    else:
        n_unloop = len(t_unloop)
    return -n_unloop * np.log(tau) - (n_cuts * t_cut.mean() + n_unloop * t_unloop.mean())/tau

def unlooping_rate_log_prior(tau, alpha=0.86, beta=9.03):
    return st.invgamma(alpha, scale=beta).logpdf(tau)
    # return st.norm(120, 2).logpdf(tau)

def unlooping_rate_log_posterior(tau, t_cut, t_unloop, **kwargs):
    log_like = unlooping_rate_log_likelihood(tau, t_cut, t_unloop)
    log_prior = unlooping_rate_log_prior(tau, **kwargs)
    return log_like  


