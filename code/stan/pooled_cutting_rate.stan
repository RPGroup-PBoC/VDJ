 
/* *****************************************************************************
* Inference of Cutting Rate Constant from Pooled Data
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* ------------------------------------------------------------------------------
* This model defines the posterior distribution for the cleavage rate constant 
* inferred from dwell times of the looped state and a determination of whether
* they lead to cutting or not. 
* ------------------------------------------------------------------------------
*/
data {
    int<lower=0> N; // Number of dwell times which lead to a cut
    int<lower=0> M; // Number of dwell times which did not cut
    vector<lower=1>[N] cut; // Dwell time in seconds for cutting events
    vector<lower=1>[M] unloop; // Dwell time in seconds for unlooping events
    }

parameters {
    real<lower=0> tau;
    }

model {
    tau ~ normal(0, 100);
    cut ~ exponential(1 / tau); 
    target += sum(-unloop / tau);
 }

generated quantities {
    real r = 1 / tau;
}
