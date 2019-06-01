/*******************************************************************************
*
* Hierarchical Model for Cutting Probability
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description
* ------------------------------------------------------------------------------
* This file defines a hierarchical model to estimate the probability of bead
* cutting for a single 12RSS mutant. The structure is as follows:
* 
*                   hyper_pcut_mu, hyper_pcut_sigma
*                                 | 
*                     day_pcut_mu, day_pcut_sigma  
*                                |
*                          replicate pcut
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
    int<lower=0> n_cuts[N]; // Number of cut beads
    int<lower=1> n_loops[N]; // Number of observed looped states
}

parameters {
    // Define the level-0 hyper parameters
    real<upper=0> log_pcut;

    // Define the level-1 hyper parameters
    vector<upper=0>[J] log_pcut_1_tilde;

    // Define the level-2 hyper parameters
    vector<upper=0>[M] log_pcut_2_tilde;

    // Define how the hyperparameters vary
    real<lower=0> tau;

}

transformed parameters { 
    // Transform parameters to non-centered versions
    vector[J] log_pcut_1 = log_pcut + tau * log_pcut_1_tilde;
    vector[M] log_pcut_2 = log_pcut_1[day_idx] + tau * log_pcut_2_tilde;

    // Transform from logspace
    real pcut = exp(log_pcut);
    vector[J] pcut_1 = exp(log_pcut_1);
    vector[M] pcut_2 = exp(log_pcut_2);
    }

model {
    // Define the hyperpriors
    log_pcut ~ normal(0, 1);
    tau ~ normal(0, 1);

    // Define the low-level priors
    log_pcut_1_tilde ~ normal(0, 1);
    log_pcut_2_tilde ~ normal(0, 1);

    // Define the likelihood
    n_cuts ~ binomial(n_loops, pcut_2[rep_idx]);
}