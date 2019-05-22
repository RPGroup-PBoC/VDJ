/*###############################################################################
# 
# Pooled Data Statistical Model 
# ------------------------------------------------------------------------------
# Author: Griffin Chure
# License: MIT
#
# Description
# ------------------------------------------------------------------------------
# This statistical model infers the five model parameters (r, r_cut, k_loop, 
# k_unloop, and sigma) ignoring day-to-day variation and pooling all data. The
# inputs are the three data sets (dwell time distribution, cutting Bernoulli
# trials, and looping fraction) as well as dimensional details. Note that this
# model has not been rigorously tested through simulation based calibration or 
# prior predictive checks, so Caveat Emptor.
*/

data {
    // Dimension specification
    int<lower=1> N; // total number of measurements

   // Data sets
   int<lower=0, upper=1> cut[N]; // Bernoulli counts of cutting
   real<lower=0, upper=1> f_looped; // Fraction of time spent in looped state
   vector<lower=0>[N] dwell_time;
}

parameters {
    real<lower=0> r; 
    real<lower=0> r_cut;
}

model {
    // Define the prior distributions
    r ~ normal(0, 0.1);
    r_cut ~ normal(0, 0.1);

    // Define the likelihoods. 
    dwell_time ~ exponential(r);
    cut ~ bernoulli(r_cut/r);
}

generated quantities {
    // Calculate k_unlooped 
    real k_unloop = r - r_cut;
    real k_loop = (r * f_looped) / (1 - f_looped);
}