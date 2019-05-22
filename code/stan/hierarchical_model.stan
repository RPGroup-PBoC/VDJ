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
    int<lower=1> J; // Number of replicates
    

    //Observables
    vector<lower=0, upper=1>[J] f_looped; // Fraction of time in looped state
    int<lower=0, upper=1> cut[N]; // Bernoulli trials for bead cutting
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

    // Homoscedastic error
    real<lower=0> sigma;
}

transformed parameters {
    // Perform non-centered transformations for better sampling
    vector<lower=0>[J] r_cut_1 = r_cut + tau * r_cut_tilde;
    vector<lower=0>[J] k_unloop_1 = k_unloop + tau * k_unloop_tilde;
    vector<lower=0>[J] k_loop_1 = k_loop + tau * k_loop_tilde;
}

model {
    // Assign the prior distributions for hyper parameters
   r_cut ~ normal(0, 1);
   k_unloop ~ normal(0, 1);
   k_loop ~ normal(0, 1);
   tau ~ normal(0, 1);

   // Prior definitions for low-level priors
   r_cut_tilde ~ normal(0, 1);
   k_unloop_tilde ~ normal(0, 1);
   k_loop_tilde ~ normal(0, 1);

   // Other priors
   sigma ~ normal(0, 0.01);

   // Evaluate the likelihoods  
   cut[idx] ~ bernoulli(r_cut_1[idx] ./ (r_cut_1[idx] + k_unloop_1[idx]));
   dwell_time[idx] ~ exponential(r_cut_1[idx] + k_unloop_1[idx]);
   f_looped ~ normal(k_loop_1 ./ (r_cut_1 + k_unloop_1 + k_loop_1), sigma);
}
