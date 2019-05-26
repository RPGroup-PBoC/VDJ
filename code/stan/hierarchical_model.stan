/*###############################################################################
#
# Hierarchical Model for Parameter Inference
# ------------------------------------------------------------------------------
# Author: Griffin Chure
# License: MIT
# 
# Description 
# ------------------------------------------------------------------------------
# This model infers the five model parameters (r, r_cut, k_loop and k_unloop)
# taking into account day-to-day variation. 
# ##############################################################################
*/

data {
    // Dimensional information
    int<lower=1> J; // Maximum number of replicates
    int<lower=1> N; // Total number of dwell time measurements
    int<lower=1> idx[N]; // ID vector for the dwell time measurements

    //Observables
    int<lower=0> total_frames[J]; // Total number of frames (per beads)
    int<lower=0> looped_frames[J]; // Total number of frames in looped states 
    int<lower=0> n_beads[J]; // Total number of beads observed per replicate 
    int<lower=0> n_cuts[J]; // Total number of cutting events per replicate
    vector<lower=0>[N] dwell_time; // dwell times for looped states
}

parameters { 
    // Hyperparameters
    real<lower=0> r_cut; 
    real<lower=0> k_unloop; 
    real<lower=0> k_loop;

    // Hyperparameters vary
    real <lower=0> tau;

    //Low level parameters
    vector<lower=0>[J] r_cut_tilde;
    vector<lower=0>[J] k_unloop_tilde;
    vector<lower=0>[J] k_loop_tilde;
}

transformed parameters {
    // Perform non-centered transformations for better sampling
    vector<lower=0>[J] r_cut_1 = r_cut + tau * r_cut_tilde;
    vector<lower=0>[J] k_unloop_1 = k_unloop + tau * k_unloop_tilde;
    vector<lower=0>[J] k_loop_1 = k_loop + tau * k_loop_tilde;
}

model {
    // Compute the cutting probabilities for each replicate; 
    vector[J] p_cut = r_cut_1 ./ (r_cut_1 + k_unloop_1);

    // Assign the prior distributions for hyper parameters
   r_cut ~ normal(0, 100);
   k_unloop ~ normal(0, 100);
   k_loop ~ normal(0, 100);
   tau ~ normal(0, 1);

   // Prior definitions for low-level priors
   r_cut_tilde ~ normal(0, 1);
   k_unloop_tilde ~ normal(0, 1);
   k_loop_tilde ~ normal(0, 1);

   // Evaluate the likelihoods  
   n_cuts ~ binomial(n_beads, p_cut);
   dwell_time[idx] ~ exponential(r_cut_1[idx] + k_unloop_1[idx]);
   looped_frames ~ binomial(total_frames, (k_loop_1 ./ (r_cut_1 + k_unloop_1 + k_loop_1)));
}
