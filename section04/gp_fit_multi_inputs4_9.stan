// Based on code from Rob Trangucci 
// https://github.com/stan-dev/example-models/tree/master/misc/gaussian-process
functions {
  matrix L_cov_exp_quad(vector[] x, 
                        real alpha,
                        vector rho,
                        real delta) {
    int N = size(x);
    matrix[N, N] K;
    real neg_half = -0.5;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(neg_half * 
                                 dot_self((x[i] - x[j]) ./ rho));
        K[j, i] = K[i, j]; 
      }
    }
    K[N, N] = sq_alpha + delta;
    return cholesky_decompose(K);
  }
}

data {
  int<lower=1> D;
  int<lower=1> N1;
  int<lower=1> N2;
  
  vector[D] X1[N1];
  vector[D] X2[N2];
  vector[N1] y;
}

transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  vector[D] X[N];
  for (n1 in 1:N1) X[n1] = X1[n1];
  for (n2 in 1:N2) X[N1 + n2] = X2[n2];
}

parameters {
  vector<lower=0>[D] rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;
}

transformed parameters {
  vector[N] f;
  {
    matrix[N, N] L_K = L_cov_exp_quad(X, alpha, rho, delta);
    f = L_K * eta; 
  }
}

model {
  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  eta ~ normal(0, 1);
  y ~ normal(f[1:N1], sigma);
}

generated quantities {
  vector[N2] y2;
  for (n2 in 1:N2)
    y2[n2] = normal_rng(f[N1 + n2], sigma);
}
