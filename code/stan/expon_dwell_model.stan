data { 
    int<lower=1> N;
    vector<lower=0>[N] dwell;
}

parameters {
    real<lower=0> tau;
}

model {
    tau ~ inv_gamma(0.86, 9.03);
    dwell ~ exponential(1/tau);
}