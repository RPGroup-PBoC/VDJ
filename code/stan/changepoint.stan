
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
    {
        vector[N+1] log_p_e;
        vector[N+1] log_p_l;
        log_p_e[1] = 0;
        log_p_l[1] = 0;
        for (i in 1:N) {
            log_p_e[i + 1] = log_p_e[i] + normal_lpdf(rms[i] | mu1, sigma1);
            log_p_l[i + 1] = log_p_l[i] + normal_lpdf(rms[i] | mu2, sigma2);
        }
    log_p = rep_vector(-log(N) + log_p_l[N+1], N) + head(log_p_e, N) -
    head(log_p_l, N);
    }
}

model {
    mu1 ~ normal(0, 100);
    mu2 ~ normal(0, 100);
    sigma1 ~ normal(0, 10);
    sigma2 ~ normal(0, 10);

    target += log_sum_exp(log_p);
}

generated quantities {
    int<lower=1, upper=N> tau; // time of changepoint
    simplex[N] sp;
    sp = softmax(log_p);
    tau = categorical_rng(sp);
}