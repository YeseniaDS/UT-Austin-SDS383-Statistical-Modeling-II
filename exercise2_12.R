# Statistical Modelling II - Section 2
rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")

#-----------------------------------------------------------
data = read.csv("dental.csv")
attach(data)
sex = as.numeric(Sex == "Male")

# add a column of 1's 
X = cbind(intercept = rep(1,nrow(data)), age, sex)
y = distance
dental = as.data.frame(cbind(age, sex, distance))
n = nrow(X); p = ncol(X)

# eyeball the plot the data 
plot(age, distance)
plot(sex, distance)

#-----------------------------------------------------------
# frequentist linear regression
# 1. least squares
lm = lm( distance ~. , data = dental)
summary(lm)
beta_lm = lm$coefficients

# 2. ridge regression
library(glmnet)
x = model.matrix(distance ~., data = dental)[,-1]
# lambdas = 10^seq(3, -2, by = -.1)
ridge = glmnet(x, y, alpha = 0, lambda = 0.1)
# cv_fit = cv.glmnet(x, y, alpha = 0, lambda = lambdas)
# plot(cv_fit)
# best_lam = cv_fit$lambda.min # cross-validation to choose best lambda
# fit = glmnet(x, y, alpha = 0, lambda = best_lam)
beta_ridge = coef(ridge)

#-----------------------------------------------------------
# Bayesian model
library(MASS)

# picking vague priors for hyperparameters
a = 2; b = 1
Lambda = diag(n)
K = diag(p) # K = diag(rep(0.01, p))
mu = matrix(0, p)

# updated parameters
Lambda_n = t(X) %*% Lambda %*% X

# for posterior beta
K_new = K + Lambda_n
mu_new = solve(K_new) %*% (t(X) %*% Lambda %*% y + K %*% mu)
# for posterior gamma
a_new = n / 2 + a
b_new = as.numeric(b + (t(y) %*% Lambda %*% y + t(mu) %*% K %*% mu 
                    - t(mu_new) %*% K_new %*% mu_new)/2)

# sampling updated parameters
nit = 10000
beta = matrix(0, nit, p)
omega = rep(0, nit)

for (i in 1:nit) {
  # omega
  tomega = rgamma(1, a_new, b_new)
  omega[i] = tomega 
  
  # beta
  beta[i, ] = mvrnorm(1, mu_new, solve(tomega * K_new))
}

# set burnin = 1000
burnin = 1001
beta_hat = colMeans(beta[burnin: nit, ])
omega_hat  = mean(omega[burnin: nit])

# Trace plots of beta's and omega
par(mar=c(4.5,4,1,3))
par(mfrow=c(2,2))
plot(omega, type = "l")
plot(beta[,1], type = "l", ylab = "beta1: intercept")
plot(beta[,2], type = "l", ylab = "beta2: age")
plot(beta[,3], type = "l", ylab = "beta3: sex = Male")

# density plot of posterior distributionn of age and sex
par(mfrow=c(1,1))
# age
hist(beta[burnin: nit,2], xlab = 'age', main = "")
abline(v=beta_lm[2], col="blue", lty=1)
abline(v=beta_ridge[2], col="red", lty=2)
abline(v=beta_hat[2], col="black", lty=3)
legend('topleft', legend=c("Least Squares", "Ridge Regression", "Bayesian"),
       col=c("blue", "red", "black"), lty=1:2, cex=0.8)

# sex = Male
hist(beta[burnin: nit,3], xlab = 'sex = Male', main = "")
abline(v=beta_lm[3], col="blue", lty=1)
abline(v=beta_ridge[3], col="red", lty=2)
abline(v=beta_hat[3], col="black", lty=3)
legend('topleft', legend=c("Least Squares", "Ridge Regression", "Bayesian"),
       col=c("blue", "red", "black"), lty=1:2, cex=0.8)
