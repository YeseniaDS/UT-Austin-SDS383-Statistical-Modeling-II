rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section05")
library(mvtnorm)
library(MASS)
library(ggplot2)

#---------------------------------------------------------------------
# Generate data points
n <- 999 
# p <- c(0.2, 0.4, 0.6, 0.8, 1)
X <- matrix(0, nrow = n, ncol = 2)
mu <- rmvnorm(5, sigma = 50 * diag(2))
count <- 0
idx <- 1
idx_sample <- matrix(0, nrow = n)
for (i in 1:n){
  # idx = sample.int(5, size = 1, prob = p)
  count <- count + 1
  if (count == 200){
    idx <- idx + 1
    count <- 0
  }
  idx_sample[i] <- idx
  X[i,] <- rmvnorm(1, mu[idx,], diag(2))
}
plot(X[,1],X[,2])

#---------------------------------------------------------------------
# Initialize
nit <- 200
c <- 10
w <- rdirichlet(1, matrix(c, nrow = 5))
mu <- mvrnorm(5, mu = colMeans(X), Sigma = diag(2))
a <- 1
b <- 1
lambda <- rgamma(n = 1, shape = a, rate = b)
mu <- as.list(data.frame(t(mu)))
D <- matrix(0, nrow = nit, ncol = nrow(X))

#---------------------------------------------------------------------
# Begin sampling
for (i in 1:nit){
  # update d
  for (j in 1:n){
    x_i <- X[j,]
    p_tmp <- mapply(function(x, w_k) w_k * 
                      dmvnorm(x_i, x, diag(2) / lambda), mu, w)
    D[i,j] <- sample(1:5, size = 1, prob = (p_tmp / sum(p_tmp)))
  }
  
  # update w
  w_param <- matrix(0, nrow = 5)
  for (k in 1:5){
    w_param[k] <- sum(D[i,] == k) + c
  }
  w <- rdirichlet(1, w_param)
  
  # update mu and lambda
  mu <- matrix(0, nrow = 5, ncol = 2)
  rate <- 0
  for (k in 1:5){
    tmp <- X[D[i,] == k, ]
    tmp <- matrix(tmp, ncol = 2)
    mu_new <- (lambda * colSums(tmp)) / (nrow(tmp) * lambda + 1)
    var_new <- diag(2) / (nrow(tmp) * lambda + 1)
    mu_tmp <- mvrnorm(n=1, mu = mu_new, Sigma = var_new)
    mu[k,] <- mu_tmp
    rate <- rate + sum((tmp - matrix(rep(mu_tmp,nrow(tmp)), ncol = 2 , byrow = T))^2)
  }
  mu <- as.list(data.frame(t(mu)))
  lambda <- rgamma(n = 1, shape = a + nrow(X), rate = b + 0.5 * rate)
  
  if (i %% 10 == 0){
    print(i)
  }
}

#---------------------------------------------------------------------
# plot: compare the original data clusters and predicted clusters
df <- data.frame(x=X[,1],y = X[,2])
ggplot(df, aes(x,y)) + geom_point(color = idx_sample) +
  ggtitle("Scatterplot of the data") + 
  theme(plot.title = element_text(size = 11))
ggplot(df, aes(x,y)) + geom_point(color = D[200,]) +
  ggtitle("Scatterplot of the labeled data by using Gibbs sampling") +
  theme(plot.title = element_text(size = 11))
