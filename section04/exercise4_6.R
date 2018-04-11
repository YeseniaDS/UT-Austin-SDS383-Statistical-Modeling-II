rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section04")
library(MASS)

#---------------------------------------------------------------------
faithful = read.csv("faithful.csv")
x <- faithful$waiting; y <- scale(faithful$eruptions)
n <- length(x)

#---------------------------------------------------------------------
# define the squared exponential function
sqexp <- function(x, alpha = 1, l = 1, sigma = 0.01) {
  
  kernel <- alpha^2 * exp(-0.5 / (l^2) * as.matrix(dist(x, upper=T, diag=T)^2))
  K.y <- kernel + sigma^2 * diag(n)
  
  return(K.y)
}

#---------------------------------------------------------------------
# define the objective function (negative log-likelihood)
loglike<- function(hypa, x, y) {
  # covareiance matrix of y
  K <- sqexp(x, alpha = hypa[1], l = hypa[2], sigma = hypa[3])
  K.inv <- solve(K)
  
  # return negative log-likelihood
  nll <- as.numeric(0.5*t(y) %*% K.inv %*% y + 
               0.5*determinant(K,logarithm=T)$modulus) + 0.5*n*log(2*pi)
  return(nll)
}

#---------------------------------------------------------------------
# initialize hyperparameters and optimize
optimal <- optim(par = c(1, 3, 1), fn = loglike, 
                 gr = NULL, method ="L-BFGS-B", x, y)
print(optimal)

#---------------------------------------------------------------------
# plot the fit
# see updated R script exercise4_5




