data { 
    int<lower=1> N; // Number of measurements
    vector<lower=0>[N] dwell; // Collection of dwell times in min
}

parameters {
    real<lower=21/60>tau1;
    real<lower=21/60, upper=tau1>tau2;
    real<lower=0, upper=1> theta;
}

model {
    tau1 ~ inv_gamma(0.86, 9.3);
    tau2 ~ inv_gamma(0.86, 9.3);
    theta ~ uniform(0, 1);
    target += sum(log((theta/tau1) * exp(-dwell ./ tau1) +
    ((1-theta)/tau2)) .* exp(-dwell ./ tau2));


}