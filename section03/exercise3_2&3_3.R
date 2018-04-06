rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
set.seed(123)
library(tidyverse)

#---------------------------------------------------------------------
# Data: survival data from the Titanic

titanic <- read_csv('titanic.csv')
summary(titanic)

titanic <- na.omit(titanic) # Age has 447 NAs

x <- matrix(titanic$Age, nrow(titanic), 1)
y <- as.numeric(titanic$Survived == "Yes")


#---------------------------------------------------------------------
# Exercise 3.2

log_posterior <- function(beta, X, y) {
  y %*% log(1/(1+exp(-X %*% beta))) + (1 - y) %*% 
    log(1 - (1/(1+exp(-X %*% beta)))) - 0.5 * t(beta) %*% beta
}

map <- optim(0, function(beta) -log_posterior(beta, x, y), 
             method = "Brent", lower = -1, upper = 1)
map$par

# [1] -0.01101471

#---------------------------------------------------------------------
# Exercise 3.3

betas <- seq(-0.2, 0.2, len = 50)
log_posterior_betas <- rep(1, 50)
for (i in 1:50) {
  log_posterior_betas[i] <- -log_posterior(betas[i], x, y)
}

par(mar=c(4.1, 4.1, 2.1, 3.1))
plot(betas, log_posterior_betas, ty= "l", xlab = expression(beta),
     ylab = "Negative Log-likelihood", col = "blue", lwd = 2, cex.main = 1,
     main = expression(paste("Posterior for Different ", beta, "s")))
points(-0.01101471, 511.4713, col = "red", pch = 19)
