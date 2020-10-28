//
// Simple mean and standard deviation computation
//

// Data
data {
  int<lower=0> N;
  vector[N] y;
  real mean_mu;
  real<lower = 0> sd_mu;
}

// Parameters
parameters {
  real mu;
  real<lower = 0> tau;
}

// Transformed parameters.
transformed parameters {
  real<lower = 0> sigma = inv_sqrt(tau);
}


// Model
model {
  mu ~ normal(mean_mu, sd_mu);
  tau ~ gamma(1e-6, 1e-6);
  y ~ normal(mu, sigma);
}

