rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
set.seed(123)

sig <- function(x) 1 / (1 + exp(-x))

log.posterior <- function(beta, X, y, lambda) {
  y %*% log(sig(X %*% beta)) + (1 - y) %*% 
    log(1 - sig(X %*% beta)) - lambda / 2 * t(beta) %*% beta
}

gradient <- function(beta, X, y, lambda) {
  t(X) %*% (y - sig(X %*% beta)) - lambda * beta
}

hessian <- function(beta, X, y, lambda) {
  - t(X) %*% diag((sig(X %*% beta) * (1 - sig(X %*% beta)))[, 1]) %*% 
    X - diag(lambda, ncol(X), ncol(X))
}

#---------------------------------------------------------------------
titanic <- read.csv("titanic.csv")
titanic <- titanic[!is.na(titanic$Age), ] # remove missing age rows

X <- as.matrix(cbind(rep(1, nrow(titanic)), scale(titanic$Age)))
y <- as.numeric(titanic$Survived == "Yes")

map <- optim(c(0, 0), function(beta)-log.posterior(beta, X, y, 1), 
             method = "L-BFGS", gr = function(beta) -gradient(beta, X, y, 1))

Mu <- map$par
Sigma <- solve(-hessian(Mu, X, y, 1))
interval <- 1.96 * sqrt(diag(Sigma))
Mu
Mu - interval
Mu + interval
