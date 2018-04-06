// Basic Poisson glm
data {
  // Define variables in data
  // Number of observations (an integer)
  int<lower=0> N;
  
  
  // Covariates
  int <lower=0, upper=1> intercept[N];
  int <lower=-1, upper=12> x[N];
  int <lower=0, upper=1> sex[N];
  
  // Count outcome
  int<lower=0> y[N];
}

parameters {
  // Define parameters to estimate
  real beta[3];
}

transformed parameters  {
  //
    real lp[N];
  real <lower=0> mu[N];
  
  for (i in 1:N) {
    // Linear predictor
    lp[i] = beta[1]*intercept[i] + beta[2]*x[i] + beta[3]*sex[i];
    
    // Mean
    mu[i] = exp(lp[i]);
  }
}

model {
  // Prior part of Bayesian inference
  beta[1]~normal(0,1);
  beta[2]~normal(0,1);
  beta[3]~normal(0,1);
  
  
  // Likelihood part of Bayesian inference
  y ~ poisson(mu);
}
