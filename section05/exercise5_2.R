rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section05")

#---------------------------------------------------------------------
data <- read.csv("restaurants.csv")
data <- data[, -1]

#---------------------------------------------------------------------
# Initialize
a <- 2
b <- 3 # prior for gamma
X <- scale(data[,1])
n <- length(X)
N <- 5000 # number of iterations

tau <- rgamma(1 ,a, 1/b)
mu_1 <- rnorm(1, 1, 1)
mu_2 <- rnorm(1, -1, 1)
mu_1_sample <- matrix(0, nrow = N)
mu_2_sample <- matrix(0, nrow = N)
sigma_sample <- matrix(1, nrow = N)
d_sample <- matrix(0, nrow = N, ncol = n)

#---------------------------------------------------------------------
# Gibbs sampling
for (i in 1:N){
  d <- matrix(1 - (runif(n) < (dnorm(X, mu_1, sqrt(1/tau)) / 
                                    (dnorm(X, mu_1, sqrt(1/tau)) +
                                       dnorm(X, mu_2, sqrt(1/tau))))), nrow = n)
  x_0 <- sum(X[d == 0])
  x_0_sq <- sum((X[d == 0])^2)
  n_0 <- sum(d == 0)
  x_1 <- sum(X[d == 1])
  x_1_sq <- sum((X[d == 1])^2)
  n_1 <- sum(d == 1)
  mu_1 <- rnorm(1, tau * x_0 / (0.1 + n_0 * tau), sqrt(1/(0.1+n_0 * tau)))
  mu_2 <- rnorm(1, tau * x_1 / (0.1 + n_1 * tau), sqrt(1/(0.1+n_1 * tau)))
  tau <- rgamma(1, a + 0.5 * (n_1 + n_0), b + (x_0_sq - 2 * x_0 * mu_1 + n_0 * mu_1 ^ 2) * 0.5 + 
                  (x_1_sq - 2 * x_1 * mu_2 + n_1 * mu_2 ^ 2) * 0.5)
  mu_1_sample[i] <- mu_1
  mu_2_sample[i] <- mu_2
  sigma_sample[i] <- sqrt(1/tau)
  d_sample[i, ] <- d
}

#---------------------------------------------------------------------
# Estimate
burnin = 1001
mu_1_est <- mean(mu_1_sample[burnin: N])
mu_2_est <- mean(mu_2_sample[burnin: N])
mean(X[data[,2]==0])
mean(X[data[,2]==1])
sigma_est <- mean(sigma_sample[burnin: N])
sum((1-d_sample[521,]) == data[,2])


# plot
hist(d_sample[1000,])
x_new <- seq(min(X), max(X), length.out = 100)
y <- 0.5 * dnorm(x_new, mu_1_est, sigma_est) + 0.5 * dnorm(x_new, mu_2_est, sigma_est)

par(mar = c(5.1, 4.1, 2.1, 2.1))
plot(density(X), main="", lwd = 2)
points(x_new, y, col = "purple", pch = 18)
legend("topright", legend = c("True density", "Estimated density"), 
       col=c("black", "purple"), lty=1:2, cex=0.7, text.font=4, bg='lightblue')

