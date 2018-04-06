setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
library(MASS)
library(ggplot2)
set.seed(51) # reproducibility

## exercise 4.3 Sampling from the prior

# define exponential kernel
kfunc_exp = function(x1, x2, l=10) {
  Sigma = matrix(rep(0, length(x1) * length(x2)), nrow = length(x1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] = exp(-0.5 * (abs(x1[i]-x2[j])/l)^2)
    }
  }
  return(Sigma)
}

# points that we will make predictions for (i.e. the test set, x.star)
x.star = matrix(seq(0, 100,len = 200), ncol = 1, nrow = 200)

# covariance matrix for the inputs
sigma = kfunc_exp(x.star, x.star, l = 0.1)

# Generate a number of functions from the process
n.samples <- 5
values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  values[,i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")

# plot the functions
ggplot(values,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="lightcyan") +
  geom_line(aes(group=variable, col=variable), lwd=1) +
  theme_bw() +
  scale_y_continuous(lim=c(-2.5,2.5), name="output, f(x)") + xlab("input, x") + 
  ggtitle("Samples from a Gaussian prior with L = 0.1") +
  theme(plot.title = element_text(hjust = 0.5))

# compare three plots of varying l = 0.1, 1, 10
l = c(0.1, 1, 10)
alpha = c(1, 5, 10)
sig1 = (alpha[1])^2 * kfunc_exp(x.star, x.star, l=10)
sig2 = (alpha[2])^2 * kfunc_exp(x.star, x.star, l=10)
sig3 = (alpha[3])^2 * kfunc_exp(x.star, x.star, l=10)
sig = matrix(list(sig1, sig2, sig3))

# Generate three functions from Gaussian process
values <- matrix(rep(0,length(x.star)*3), ncol=3)
for (i in 1:3) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  sigma = sig[[i]]
  values[,i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")

ggplot(values,aes(x=x,y=value)) +
  geom_line(aes(group=variable, col=variable), lwd=1) +
  theme_bw() +
  scale_y_continuous(lim=c(-14,14), name="output, f(x)") + xlab("input, x") + 
  ggtitle("Samples from a Gaussian prior with Varying Alpha") +
  theme(plot.title = element_text(hjust = 0.5))

