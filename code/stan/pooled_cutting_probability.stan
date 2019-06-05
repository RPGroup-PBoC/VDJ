/*******************************************************************************
*
* Pooled Model for Cutting Probability
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description
* ------------------------------------------------------------------------------
* This file defines a single-level model to estimate the probability of bead
* cutting for a single 12RSS mutant. 
* ------------------------------------------------------------------------------
*/
data {
    int<lower=0> n_cuts; // Number of cut beads
    int<lower=1> n_loops; // Number of observed looped states
}

parameters {
    // Define the level-0 hyper parameters
    real<upper=0> log_pcut;
}

transformed parameters { 
    // Transform from logspace
    real pcut = exp(log_pcut);
    }

model {
    // Define the hyperpriors
    log_pcut ~ normal(0, 4);
    // Define the likelihood
    n_cuts ~ binomial(n_loops, pcut);
}