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
    int<lower=1> J; // Number of replicates, used for portability of helper code 

   // Data sets
   int<lower=0> n_beads[J]; // Total number of beads 
   int<lower=0> n_cuts[J]; // Number of cut beads
   int<lower=0> total_frames[J]; // Total number of frames
   int<lower=0> looped_frames[J];
   vector<lower=0>[N] dwell_time;
}

transformed data {
    // Remove replicate information
    int<lower=0> n_beads_reps = sum(n_beads);
    int<lower=0> n_cut_reps = sum(n_cuts);
    int<lower=0> total_frames_reps = sum(total_frames);
    int<lower=0> looped_frames_reps = sum(looped_frames); 
}

parameters {
    real<lower=0>  r_cut;
    real<lower=0> k_loop;
    real<lower=0> k_unloop;
}

model {
    // Define the prior distributions
    r_cut ~ normal(0, 100);
    k_loop ~ normal(0, 100);
    k_unloop ~ normal(0, 100);

    // Define the likelihoods. 
    dwell_time ~ exponential(r_cut + k_unloop);
    n_cut_reps ~ binomial(n_beads_reps, (r_cut / (r_cut + k_unloop))); 
    looped_frames_reps ~ binomial(total_frames_reps, (k_loop / (r_cut + k_loop + k_unloop)));
}
