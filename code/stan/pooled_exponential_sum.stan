// functions {
//     real sum_expon_lpdf(vector time, real theta, positive_ordered tau) {
//         return sum(log(theta/tau[1]) - time./tau[1] + log((1-theta)/tau[2]) - time./tau[2]);
//     }
// }
data { 
    int<lower=1> N; // Number of measurements
    vector<lower=0>[N] dwell; // Collection of dwell times in min
}

parameters {
    positive_ordered[2] tau;
    real<lower=0> theta;
}

model {
    tau ~ inv_gamma(0.86, 9.3);
    theta ~ normal(0, 0.1);
    target += sum(log((theta/tau[1]) * exp(-dwell ./ tau[1]) +
    ((1-theta)/tau[2])) .* exp(-dwell ./ tau[2]));


}