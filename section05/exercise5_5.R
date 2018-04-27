rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section05")

#---------------------------------------------------------------------
data = read.csv("restaurants.csv")
data = data[, -1]
x = scale(data[,1])
n = length(x)

#---------------------------------------------------------------------
# Parameters
a = 1; b = 1  # inverse gamma parameters
m_0 = 1; m_1 = -1 # mean of normal prior for mu0 and mu1
l = 1/100 # prior precision for means

# Initialize
N = 10000 # number of iterations
xnew = numeric(0)
mu0 = numeric(0); mu1 = numeric(0)
mu0[1] = -0.1; mu1[1] = 0.1
lambda = numeric(0); lambda[1] = 1 # precision of the model
sigma2 = numeric(0); sigma2[1] = 1/lambda[1]
w = numeric(0)
n0 = numeric(0); n1 = numeric(0)
d = matrix(0, n, N+1)
r = numeric(0)
R = matrix(0, n, N)

# set x values to evaluate density
x0 = min(x); x1 = max(x)
d0 = x0 - 0.1 * (x1-x0); d1 = x1 + 0.1 * (x1-x0)
dx = seq(d0, d1, length = 10 * n)
dens.matrix = matrix(0, ncol = length(dx), nrow = N)

#---------------------------------------------------------------------
# Gibbs sampling
for(k in 1:N) {
  n1[k] = sum(d[,k])
  n0[k] = n - n1[k]
  
  w[k+1] = rbeta(1, a + n0[k], b + n1[k])
  
  m0 = (l * m_0 * sigma2[k] + sum((1 - d[,k]) * x)) / (l*sigma2[k] + n0[k])
  m1 = (l * m_1 * sigma2[k] + sum(d[,k] * x))/ (l*sigma2[k] + n1[k])
  mu0[k+1] = rnorm(1, m0, (l + n0[k]/sigma2[k])^(-1/2))
  mu1[k+1] = rnorm(1, m1, (l + n1[k]/sigma2[k])^(-1/2))
  
  for (i in 1:n) {
    r[i] = w[k+1] * dnorm(x[i], mu0[k+1], sqrt(sigma2[k]))/ 
      (w[k+1] * dnorm(x[i], mu0[k+1], sqrt(sigma2[k])) + 
         (1-w[k+1]) * dnorm(x[i], mu1[k+1], sqrt(sigma2[k])))
    d[i, k+1] = 1 - rbinom(1, 1, r[i])
    R[i,k] = r[i]
  }
  
  # component precision
  lambda[k+1] = rgamma(1, a + n/2, b + (sum((1 - d[,k+1]) * ((x-mu0[k+1])^2)) 
                                        + sum(d[,k+1] * ((x-mu1[k+1])^2)))/2)
  sigma2[k+1] = 1/lambda[k+1]
  
  # sample x
  y = w[k+1] * dnorm (dx, mu0[k+1], sqrt(sigma2[k+1])) +
    (1-w[k+1]) * dnorm (dx, mu1[k+1], sqrt(sigma2[k+1])) 
  dens.matrix[k, ] = y
}

#---------------------------------------------------------------------
burnin = 2001
mu0_est <- mean(mu0[burnin: N])
mu1_est <- mean(mu1[burnin: N])
mean(x[data[,2]==0])
mean(x[data[,2]==1])
sigma2_est <- mean(sigma2[burnin: N])
weight_est <- mean(w[burnin: N])

# plots
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 2.1, 3.1))
# trace plots of mu0, mu1, sigma2, and w
plot(mu0[burnin:N], type="l", xlab = "", ylab = expression(mu[0]), 
     main = "Gibbs Sampler Simulation")
plot(mu1[burnin:N], type="l", xlab="Number of Iterations", 
     ylab = expression(mu[1]), main = "")
plot(sigma2[burnin:N], type="l", xlab = "", ylab = expression(sigma^2), 
     main = "Gibbs Sampler Simulation")
plot(w[burnin:N], type="l", xlab = "Number of Iterations", ylab = "w")

# histogram of mu0 and mu1, sigma2, and weights w
par(mfrow = c(2, 2))
hist(mu0[burnin:N], xlab = expression(mu[0]), cex.main = 1,
     main = expression(paste("Histogram  of ", mu [0])))
abline(v = mean(mu0[burnin:N]), col="blue")
hist(mu1[burnin:N], xlab = expression(mu[1]), cex.main = 1,
     main = expression(paste("Histogram  of ", mu [1])))
abline(v = mean(mu1[burnin:N]), col="blue")
hist(sigma2[burnin:N], xlab = expression(sigma^2),cex.main = 1,
     main = expression(paste("Histogram  of ", sigma^2)))
abline(v= sigma2_est, col="blue")
hist(w[burnin:N], xlab = "w", cex.main = 1, 
     main = expression(paste("Histogram  of ", pi)))
abline(v= mean(w[burnin:N]), col="blue")

# compare the original density with the estimated density
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 2.1, 3.1))
g.hat = colMeans(dens.matrix) 
hist(x, freq = F, breaks = 25, main = "")
lines(density(x), main="", lwd = 2)
lines(dx, g.hat, col = "purple", lwd = 3)
legend("topright", legend = c("True density", "Estimated density"), 
       col=c("black", "purple"), lty=1:2, cex=0.7, text.font=4, bg='lightblue')

