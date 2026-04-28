// Two-sample AR(1) hierarchical model for htdv_simstudy().
// Shared autocorrelation and innovation scale; group-specific means.
// Test parameter: delta = alpha2 - alpha1.

data {
  int<lower=1> N1;
  int<lower=1> N2;
  vector[N1] x1;
  vector[N2] x2;
}
parameters {
  real alpha1;
  real alpha2;
  real<lower=-0.999, upper=0.999> phi;
  real<lower=0> sigma;
}
transformed parameters {
  real delta;
  delta = alpha2 - alpha1;
}
model {
  alpha1 ~ normal(0, 10);
  alpha2 ~ normal(0, 10);
  phi ~ normal(0, 0.5);
  sigma ~ student_t(3, 0, 2.5);

  x1[1] ~ normal(alpha1, sigma / sqrt(1 - square(phi)));
  for (t in 2:N1)
    x1[t] ~ normal(alpha1 + phi * (x1[t-1] - alpha1), sigma);

  x2[1] ~ normal(alpha2, sigma / sqrt(1 - square(phi)));
  for (t in 2:N2)
    x2[t] ~ normal(alpha2 + phi * (x2[t-1] - alpha2), sigma);
}
