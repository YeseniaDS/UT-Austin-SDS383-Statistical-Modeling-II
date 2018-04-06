rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
set.seed(123)
library(dplyr)

#---------------------------------------------------------------------
# Data: diabetes information for women of Pima indian heritage

pima <- read.csv('pima.csv')
summary(pima)
# pairs(pima)

y <- as.matrix(pima$class_variable )
n <- nrow(y)

X <- as.matrix(cbind(intercept = rep(1,n), pima[,-9]))
p <- ncol(X)

#---------------------------------------------------------------------
# Frequentist Estimation

# Fit the model to the data
fit <- glm(y ~ X, family = binomial(link = probit))
beta.mle <- fit$coefficients
print(beta.mle)

#---------------------------------------------------------------------
# Bayesian Estimation - Gibb Sampler

library(MASS)
library(truncnorm) # Truncated Normal distribution

# Prior Hyperparameters
K.zero <- diag(rep(0.1, p)) # precision matriz pxp
beta.zero <- matrix(0, p) # prior guess on beta
a.zero <- 1 # prior sample size for the error variance
b.zero <- 1 # prior sum of square errors for the error variance

N1  <- sum(y)  # Number of successes
N0  <- n - N1  # Number of failures

# Sampling updated parameters
n.iter <- 5000
beta <- matrix(NA, n.iter, p)
tau <- matrix(NA, n.iter)
beta[1,] <- rep(0, p)
tau[1] <- 1
z <- rep(0, n)

# Gibb Sampler
for (i in 2:n.iter) {
  # Update mean of z based on beta
  mu_z <- X %*% beta[i-1,]
  # Get the latent variable from its distribution
  z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
  z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
  
  # Get the Betas
  K.new <- K.zero + crossprod(X)
  beta.new <- solve(K.new) %*% (crossprod(X, z) + K.zero %*% beta.zero)
  beta[i,] <- mvrnorm(1, beta.new, solve(tau[i-1] * K.new))
  
  # Get the tau
  a.new <- a.zero + (n + 1) / 2
  s <- t(beta.zero) %*% K.zero %*% beta.zero + crossprod(z)
  r <- t(beta.zero) %*% (K.zero + crossprod(X)) %*% beta.zero
  b.new <- as.numeric(b.zero + 1 / 2 * (s - r))
  tau[i] <- rgamma(1, a.new, b.new)
}

burn.in <- 1000
beta.post <- colMeans(beta[-(1:burn.in), ])
y.hat <- ifelse(X %*% beta.post > 0, 1, 0)
accuracy.bayes <- sum(y.hat)/sum(y) * 100

# Trace Plots
par(mar=c(4.1, 4.1, 2.1, 3.1))
par(mfrow=c(4,2))
for (i in 1:8) {
  plot(beta[1001:5000,i], xlab = "Iterations", 
       ylab = expression(beta[i]), type = "l")
}
