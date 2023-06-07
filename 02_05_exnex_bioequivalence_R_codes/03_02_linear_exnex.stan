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
  real m_mean_mu;
  real<lower = 0> t_mean_mu;
  real<lower = 0> shape_sigma_mu;
  real<lower = 0> scale_sigma_mu;
}

// Parameter
parameters {
  real mu_tilde[N_formulation];
  real mean_mu;
  real<lower = 0> sigma_mu;
  real beta;
  real<lower = 0> sigma_epsilon;
  real<lower = 0> sigma_s;
  real s[N_subjid];
  real<lower = 0, upper = 1> pi[N_formulation];
}

transformed parameters {
  real mu[N_formulation];
  real mu_y[N];

  for (i in 1:N_formulation) {
    mu[i] = mean_mu + sigma_mu * mu_tilde[i];
  }
  for (n in 1:N) {
    mu_y[n] = mu[index_formulation[n]] + beta * period_covariate[n] + s[index_subjid[n]];
  }
}

// Model
model {
  sigma_epsilon ~ normal(0.0, 100.0);
  sigma_s ~ normal(0.0, 100.0);
  mu_tilde ~ normal(0.0, 1.0);
  mean_mu ~ normal(m_mean_mu, t_mean_mu);
  sigma_mu ~ inv_gamma(shape_sigma_mu, scale_sigma_mu);
  beta ~ normal(0.0, 100.0);
  s ~ normal(0.0, sigma_s);
  for (n in 1:N){
    target += log_mix(pi[],
                      normal_lpdf(y[n]|theta_ex[index_formulation[n]], sigma[]),
                      normal_lpdf(y[n]|theta_nonex[index_formulation[n]], sigma[]));
  }
}
