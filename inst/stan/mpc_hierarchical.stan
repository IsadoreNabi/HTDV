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
  real<lower=0> psi0;
}
model {
  target += normal_lpdf(theta | mu0, tau0);
  target += half_cauchy_lpdf(sigma_inf | sigma_inf_scale);
  target += gamma_lpdf(gamma_mix - 1 | gamma_shape, gamma_rate);
  target += half_cauchy_lpdf(psi0 | sigma_inf_scale);
  {
    real var_tot = psi0 * psi0 *
      (gamma_mix / (gamma_mix - 1))^2 * (1 + 4 / (gamma_mix - 1));
    real scale_step = sqrt(var_tot / N);
    target += normal_lpdf(x[1] | theta, scale_step);
    for (t in 2:N) {
      real weight = pow(t + 0.0, -gamma_mix / 2.0);
      target += normal_lpdf(x[t] | theta + weight * (x[t - 1] - theta),
                            scale_step);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] x_rep;
  {
    real var_tot = psi0 * psi0 *
      (gamma_mix / (gamma_mix - 1))^2 * (1 + 4 / (gamma_mix - 1));
    real scale_step = sqrt(var_tot / N);
    log_lik[1] = normal_lpdf(x[1] | theta, scale_step);
    x_rep[1] = normal_rng(theta, scale_step);
    for (t in 2:N) {
      real weight = pow(t + 0.0, -gamma_mix / 2.0);
      log_lik[t] = normal_lpdf(x[t] | theta + weight * (x[t - 1] - theta),
                               scale_step);
      x_rep[t] = normal_rng(theta + weight * (x_rep[t - 1] - theta),
                            scale_step);
    }
  }
}
