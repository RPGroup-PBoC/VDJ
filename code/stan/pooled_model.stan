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
    int<lower=1> N; // total number of dwell times;

   // Data sets
   int<lower=0> n_beads; // Total number of beads 
   int<lower=0> n_cuts; // Number of cut beads
   int<lower=0> total_frames; // Total number of frames
   int<lower=0> looped_frames;
   vector<lower=0>[N] dwell_time;
}

parameters {
    real<lower=0>  r_cut;
    real<lower=0> k_loop;
    real<lower=0> k_unloop;
}

model {
    // Define the prior distributions
    r_cut ~ normal(0, 1);
    k_loop ~ normal(0, 1);
    k_unloop ~ normal(0, 1);

    // Define the likelihoods. 
    dwell_time ~ exponential(r_cut + k_unloop);
    n_cuts ~ binomial(n_beads, (r_cut / (r_cut + k_unloop))); 
    looped_frames ~ binomial(total_frames, (k_loop / (r_cut + k_loop + k_unloop)));
}

generated quantities {
    real p_loop = k_loop / (k_loop + k_unloop);
    real p_cut = r_cut / (r_cut + k_unloop);
    real f_loop = k_loop / (k_loop + k_unloop + r_cut);
}