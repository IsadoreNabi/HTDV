functions {
  real half_cauchy_lpdf(real y, real scale) {
    return cauchy_lpdf(y | 0, scale) + log(2);
  }
  real pairwise_gauss_lpdf(vector pair, real theta, real sigma_inf, real phi) {
    real a = pair[1];
    real b = pair[2];
    real sigma_marg = sigma_inf;
    real sigma_cond = sigma_inf * sqrt(1 - phi * phi);
    return normal_lpdf(a | theta, sigma_marg)
           + normal_lpdf(b | theta + phi * (a - theta), sigma_cond);
  }
}
data {
  int<lower=3> N;
  vector[N] x;
  real<lower=0> sigma_inf_scale;
  real<lower=0> gamma_shape;
  real<lower=0> gamma_rate;
  real mu0;
  real<lower=0> tau0;
}
parameters {
  real theta;
  real<lower=0> sigma_inf;
  real<lower=1> gamma_mix;
  real<lower=-0.999, upper=0.999> phi;
}
model {
  target += normal_lpdf(theta | mu0, tau0);
  target += half_cauchy_lpdf(sigma_inf | sigma_inf_scale);
  target += gamma_lpdf(gamma_mix - 1 | gamma_shape, gamma_rate);
  for (t in 1:(N - 1)) {
    vector[2] pair;
    pair[1] = x[t];
    pair[2] = x[t + 1];
    target += pairwise_gauss_lpdf(pair | theta, sigma_inf, phi);
  }
}
generated quantities {
  vector[N - 1] log_lik;
  vector[N] x_rep;
  real sigma_cond = sigma_inf * sqrt(1 - phi * phi);
  for (t in 1:(N - 1)) {
    vector[2] pair;
    pair[1] = x[t];
    pair[2] = x[t + 1];
    log_lik[t] = pairwise_gauss_lpdf(pair | theta, sigma_inf, phi);
  }
  x_rep[1] = normal_rng(theta, sigma_inf);
  for (t in 2:N) {
    x_rep[t] = normal_rng(theta + phi * (x_rep[t - 1] - theta), sigma_cond);
  }
}
