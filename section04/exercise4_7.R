rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section04")
library(MASS)

#---------------------------------------------------------------------
faithful = read.csv("faithful.csv")
# index <- sample(1:nrow(faithful), 10)
index <- c(105, 68, 97, 72, 91, 93, 59, 156, 246, 139)
x <- faithful$waiting[index]
y <- scale(faithful$eruptions[index])
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
alpha <- c((1:10)/10, 2:100)
l <- c((1:10)/10, 2:100)
sigma <- c((1:10)/100, (2:10)/10, 2:91)

niter <- 80
parmat <- matrix(NA, nrow = niter, ncol = 3)
for (i in 1:niter){
  par <- c(sample(alpha, 1), sample(l, 1), sample(sigma, 1))
  optimal <- optim(par = par, fn = loglike, lower = 0.01,
                   gr = NULL, method ="L-BFGS-B", x, y)
  parmat[i, ] <- optimal$par
  
}

#---------------------------------------------------------------------
# plot the histogram of hyper-parameters
par(mfrow = c(3, 2))
par(mar = c(4.1, 4.0, 1.1, 2.1))
plot(parmat[,1], ylab = expression(alpha), pch = 19, cex = 0.8, col = "red")
plot(parmat[,2], ylab = "length-scale l", pch = 19, cex = 0.8, col ="blue")
plot(parmat[,3], ylab = expression(sigma), pch = 19, cex = 0.8, col ="orange")
hist(parmat[,2], xlab = "length-scale l", main = "", breaks = 20)
hist(parmat[,1], xlab = expression(alpha), main = "", breaks = 20)
hist(parmat[,3], xlab = expression(sigma), main = "", breaks = 20)








