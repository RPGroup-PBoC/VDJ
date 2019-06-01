/*******************************************************************************
*
* Hierarchical Model for Fraction of Time Looped
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description
* ------------------------------------------------------------------------------
* This file defines a hierarchical model to estimate the probability of being in
* the looped state for a single 12RSS mutant.
* ------------------------------------------------------------------------------
*/
data {
    // Dimensional parameters 
    int<lower=1> J; // Total number of days
    int<lower=1> M; // Number of unique replicates
    int<lower=1> N; // Total number of measurements

    // Hierarchy identification vectors
    int<lower=1, upper=J> day_idx[M]; // Identifies the day 
    int<lower=1, upper=M> rep_idx[N]; // Identifies the replicate

    // Observed parameters
    int<lower=0> total_frames[N]; // Number of cut beads
    int<lower=0> looped_frames[N]; // Number of observed looped states
}

parameters {
    // Define the level-0 hyper parameters
    real<upper=0> log_ploop;

    // Define the level-1 hyper parameters
    vector<upper=0>[J] log_ploop_1_tilde;

    // Define the level-2 hyper parameters
    vector<upper=0>[M] log_ploop_2_tilde;

    // Define how the hyperparameters vary
    real<lower=0> tau;

}

transformed parameters { 
    // Transform parameters to non-centered versions
    vector[J] log_ploop_1 = log_ploop + tau * log_ploop_1_tilde;
    vector[M] log_ploop_2 = log_ploop_1[day_idx] + tau * log_ploop_2_tilde;

    // Transform from logspace
    real ploop = exp(log_ploop);
    vector[J] ploop_1 = exp(log_ploop_1);
    vector[M] ploop_2 = exp(log_ploop_2);
    }

model {
    // Define the hyperpriors
    log_ploop ~ normal(0, 1);
    tau ~ normal(0, 1);

    // Define the low-level priors
    log_ploop_1_tilde ~ normal(0, 0.1);
    log_ploop_2_tilde ~ normal(0, 0.1);

    // Define the likelihood
    looped_frames ~ binomial(total_frames, ploop_2[rep_idx]);
}