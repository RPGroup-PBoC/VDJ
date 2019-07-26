
data {
    int<lower=1> N; // Number of time measurements
    real rms[N]; // Root mean squared displacement
}

parameters {
    real mu1;
    real mu2;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
}

transformed parameters { 
    vector[N] log_p;
    real mu;
    real sigma;
    log_p = rep_vector(-log(N), N);
    for (t in 1:N) {
        for (i in 1:N) {
            mu = i < t ? mu1 : mu2;
            sigma = i < t ? sigma1: sigma2;
            log_p[t] += normal_lpdf(rms[i] | mu, sigma);
        }
    }
}

model {
    mu1 ~ normal(0, 100);
    mu2 ~ normal(0, 100);
    sigma1 ~ normal(0, 100);
    sigma2 ~ normal(0, 100);

    target += log_sum_exp(log_p);
}

generated quantities {
    int<lower=1, upper=N> tau; // time of changepoint
    simplex[N] sp;
    sp = softmax(log_p);
    tau = categorical_rng(sp);
}