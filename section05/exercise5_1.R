rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section05")

#---------------------------------------------------------------------
dat = read.csv("restaurants.csv", sep = ",", header = TRUE)[, -1]

smp_size = floor(0.75 * nrow(dat))
set.seed(1) # reproducibility
train_ind = sample(seq_len(nrow(dat)), size = smp_size)

# train dataset and test dataset
train = dat[train_ind, ]
test = dat[-train_ind, ]

#---------------------------------------------------------------------
# fit a linear model
lmfit = lm(Profit ~ factor(DinnerService) + SeatingCapacity, data = train)
print(lmfit$coefficients)

y.hat = predict(lmfit, test)
rmse = sqrt(mean((test$Profit-y.hat)^2))
print(rmse)

#---------------------------------------------------------------------
# Bayesian model
library(MASS)
n = nrow(train)
p = ncol(train)
X = as.matrix(cbind(intercept = rep(1, n), train[,2:3]))
y = train$Profit

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

# compute the rmse
error <- sqrt(as.numeric(crossprod(y - X %*% beta_hat)/(n-p)))
print(error)

#---------------------------------------------------------------------
# Traceplots of beta's and omega
par(mar=c(4.5,4,1,3))
par(mfrow=c(2,2))
plot(omega, type = "l")
plot(beta[,1], type = "l", ylab = "Intercept")
plot(beta[,2], type = "l", ylab = "Dinner Service ")
plot(beta[,3], type = "l", ylab = "Seating Capacity")

# plot the residulas
par(mfrow=c(1,2))
err <- y - X %*% beta_hat
hist(err, xlab = "Residuals", freq = F, breaks = 25,
     main = "Histogram of Residuals", cex.main = 1)
hist(train$Profit, xlab = "Profit", freq = F, breaks = 25,
     main = "Histogram of Profit", cex.main = 1)

