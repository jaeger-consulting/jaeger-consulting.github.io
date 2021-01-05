//
// Linear 2x2 crossover model
//

// Data
data {
  int<lower = 0> N;
  int<lower = 0> N_formulation;
  int<lower = 0> N_subjid;
  real y[N];
  int<lower = 1, upper = N_formulation> index_formulation[N];
  real period_covariate[N];
  int<lower = 1, upper = N_subjid> index_subjid[N];
}

// Parameter
parameters {
  real mu[N_formulation];
  real beta;
  real<lower = 0> sigma_epsilon;
  real<lower = 0> sigma_s;
  real s[N_subjid];
}

transformed parameters {
  real mu_y[N];

  for (n in 1:N) {
    mu_y[n] = mu[index_formulation[n]] + beta * period_covariate[n] + s[index_subjid[n]];
  }
}

// Model
model {
  sigma_epsilon ~ normal(0.0, 100.0);
  sigma_s ~ normal(0.0, 100.0);
  mu ~ normal(0.0, 100.0);
  beta ~ normal(0.0, 100.0);
  s ~ normal(0.0, sigma_s);
  y ~ normal(mu_y, sigma_epsilon);
}
