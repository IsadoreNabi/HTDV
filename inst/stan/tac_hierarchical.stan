functions {
  real half_cauchy_lpdf(real y, real scale) {
    return cauchy_lpdf(y | 0, scale) + log(2);
  }
}
data {
  int<lower=2> N;
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
transformed parameters {
  real<lower=0> sigma_e = sigma_inf * sqrt(1 - phi * phi);
}
model {
  target += normal_lpdf(theta | mu0, tau0);
  target += half_cauchy_lpdf(sigma_inf | sigma_inf_scale);
  target += gamma_lpdf(gamma_mix - 1 | gamma_shape, gamma_rate);
  target += normal_lpdf(x[1] | theta, sigma_inf);
  for (t in 2:N) {
    target += normal_lpdf(x[t] | theta + phi * (x[t - 1] - theta), sigma_e);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] x_rep;
  log_lik[1] = normal_lpdf(x[1] | theta, sigma_inf);
  x_rep[1] = normal_rng(theta, sigma_inf);
  for (t in 2:N) {
    log_lik[t] = normal_lpdf(x[t] | theta + phi * (x[t - 1] - theta), sigma_e);
    x_rep[t] = normal_rng(theta + phi * (x_rep[t - 1] - theta), sigma_e);
  }
}
