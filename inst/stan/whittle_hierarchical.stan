functions {
  real half_cauchy_lpdf(real y, real scale) {
    return cauchy_lpdf(y | 0, scale) + log(2);
  }
  real ar1_spectrum(real lambda, real phi, real sigma_e) {
    real num = sigma_e * sigma_e;
    real den = 2 * pi() * (1 - 2 * phi * cos(lambda) + phi * phi);
    return num / den;
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
transformed data {
  int M = N / 2;
  vector[M] lambda_grid;
  vector[M] periodogram;
  real xbar = mean(x);
  {
    vector[N] xc = x - xbar;
    vector[N] re = rep_vector(0.0, N);
    vector[N] im = rep_vector(0.0, N);
    for (j in 1:M) {
      real lam = 2 * pi() * j / N;
      lambda_grid[j] = lam;
      real re_sum = 0;
      real im_sum = 0;
      for (t in 1:N) {
        re_sum += xc[t] * cos(lam * (t - 1));
        im_sum -= xc[t] * sin(lam * (t - 1));
      }
      periodogram[j] = (re_sum * re_sum + im_sum * im_sum) / (2 * pi() * N);
    }
  }
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
  target += normal_lpdf(xbar | theta, sigma_inf / sqrt(N + 0.0));
  for (j in 1:M) {
    real f = ar1_spectrum(lambda_grid[j], phi, sigma_e);
    target += -log(f) - periodogram[j] / f;
  }
}
generated quantities {
  vector[M] log_lik;
  vector[N] x_rep;
  for (j in 1:M) {
    real f = ar1_spectrum(lambda_grid[j], phi, sigma_e);
    log_lik[j] = -log(f) - periodogram[j] / f;
  }
  x_rep[1] = normal_rng(theta, sigma_inf);
  for (t in 2:N) {
    x_rep[t] = normal_rng(theta + phi * (x_rep[t - 1] - theta), sigma_e);
  }
}
