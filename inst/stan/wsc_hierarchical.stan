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
  int<lower=1> b_lower;
  int<lower=1> b_upper;
}
parameters {
  real theta;
  real<lower=0> sigma_inf;
  real<lower=1> gamma_mix;
  real<lower=-0.999, upper=0.999> phi;
  real<lower=0> log_b_raw;
}
transformed parameters {
  real<lower=0> sigma_e = sigma_inf * sqrt(1 - phi * phi);
  real<lower=b_lower, upper=b_upper> b_cont;
  b_cont = b_lower +
    (b_upper - b_lower) * inv_logit(log(log_b_raw + 1e-6) - log(1 - log_b_raw + 1e-6));
}
model {
  target += normal_lpdf(theta | mu0, tau0);
  target += half_cauchy_lpdf(sigma_inf | sigma_inf_scale);
  target += gamma_lpdf(gamma_mix - 1 | gamma_shape, gamma_rate);
  target += beta_lpdf(log_b_raw | 2, 2);
  target += normal_lpdf(x[1] | theta, sigma_inf);
  for (t in 2:N) {
    real local_phi = phi * exp(-fmax(0.0, (1.0 / b_cont - 0.01)));
    target += normal_lpdf(x[t] | theta + local_phi * (x[t - 1] - theta), sigma_e);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] x_rep;
  log_lik[1] = normal_lpdf(x[1] | theta, sigma_inf);
  x_rep[1] = normal_rng(theta, sigma_inf);
  {
    real local_phi = phi * exp(-fmax(0.0, (1.0 / b_cont - 0.01)));
    for (t in 2:N) {
      log_lik[t] = normal_lpdf(x[t] | theta + local_phi * (x[t - 1] - theta),
                               sigma_e);
      x_rep[t] = normal_rng(theta + local_phi * (x_rep[t - 1] - theta),
                            sigma_e);
    }
  }
}
