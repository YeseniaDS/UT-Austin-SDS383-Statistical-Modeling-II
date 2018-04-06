rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
set.seed(123)

#---------------------------------------------------------------------
log.posterior <- function(beta, X, y, lambda) {
  t(y) %*% X %*% beta - sum(exp(X %*% beta)) - lambda / 2 * t(beta) %*% beta
}

gradient <- function(beta, X, y, lambda) {
  as.vector(t(X) %*% (y - exp(X %*% beta)) - lambda * beta)
}

hessian <- function(beta, X, y, lambda) {
  -t(X) %*% diag(as.vector(exp(X %*% beta))) %*% X - lambda
}


#---------------------------------------------------------------------
tea <- read.csv("tea_discipline_oss.csv", stringsAsFactors = F)
tea <- tea[tea$ACTIONS > 0, ]

X <- as.matrix(cbind(rep(1, nrow(tea)), scale(as.numeric(tea$GRADE))))
y <- tea$ACTIONS

map <- optim(c(0, 0), function(x) -log.posterior(x, X, y, 1), 
             gr = function(x) -gradient(x, X, y, 1))
# check: glm(y ~ X -1, family = poisson(link = "log"))

Mu <- map$par
Cov <- solve(-hessian(Mu, X, y, 1))
interval <- 1.96 * sqrt(diag(Cov))
Mu 
Mu - interval
Mu + interval