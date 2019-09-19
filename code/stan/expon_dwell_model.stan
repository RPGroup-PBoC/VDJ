/* *****************************************************************************
* Exponential Dwell Time Model
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* ------------------------------------------------------------------------------
* This model infers the shape parameter tau of and exponential distribution to 
* dwell time meausurements of the paired complex. The prior on tau is defined as
* a beta distribution with support on time in units of minutes. The provided 
* dwell time measurements should therefore be provided in minutes. 
* *****************************************************************************/

data { 
    int<lower=1> N; // Number of measurements
    vector<lower=0>[N] dwell; // Dwell time measurements in units of min
}

parameters {
    real<lower=0> tau; // Average time between arrivals of paired complex dissociation
}

model {
    // Prior determined bye eye
    tau ~ inv_gamma(0.86, 9.03);

    // Likelihood 
    dwell ~ exponential(1/tau);
}