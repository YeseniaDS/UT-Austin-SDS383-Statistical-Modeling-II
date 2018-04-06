setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
library(MASS)
library(ggplot2)
set.seed(51) # reproducibility

## exercise 4.5 

# define exponential kernel
kfunc_exp = function(x1, x2, l=1) {
  Sigma = matrix(rep(0, length(x1) * length(x2)), nrow = length(x1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] = exp(-0.5 * (abs(x1[i]-x2[j])/l)^2)
    }
  }
  return(Sigma)
}

# function to get the GP solution
gp_solve = function(x.train, y.train, x.pred, kernel, sigma2e = 0) {
  solution = list()
  # compute training covariance matrix (used to get relationships in training)
  k.xx = kernel(x.train,x.train)
  # compute covariance between training and testing (used to predict weights into new data)
  k.x_xp = kernel(x.train,x.pred)
  k.xp_x = kernel(x.pred,x.train)
  # compute covariance between testing (used to estimate covariance of predictions in new data)
  k.xp_xp = kernel(x.pred,x.pred)
  
  # compute the prediction without noise = K(xs,x) K(x,x)^{-1} y [eq 2.19]
  # Vinv = solve(k.xx)
  # Mean and covariance functions now include noise [eq 2.22-2.24], 
  # which only requires adding an identity term to the variance component
  Vinv = ginv(k.xx + sigma2e * diag(1, ncol(k.xx)))
  
  # compute the estimate and variance for the prediction [eq 2.19 or 2.23]
  solution[["mu"]] = k.xp_x %*% Vinv %*% y.train
  solution[["var"]] = k.xp_xp - k.xp_x %*% Vinv %*% k.x_xp
  return( solution )
}

# function to plot the GP solution
gp_plot = function(x, y, x.pred, mu.pred, cov.pred, n.sample=0, main="", scale.se=T){
  if (scale.se) {
    # scale around credible intervals
    yrng = range(c(y, mu.pred + 1.96 * sqrt(diag(cov.pred)), 
                   mu.pred - 1.96 * sqrt(diag(cov.pred))), na.rm = T)
  } else {
    # scale around training data (useful for multiple plots of same data)
    yrng = range(y, na.rm = T) * 1.05
  }
  
  xrng = range(c(x,x.pred), na.rm = T)
  
  plot(x, y, type="n", xlab="x", ylab="y", 
       xlim=xrng, ylim=yrng, xaxs="i", main=main)
  
  # plot the credible interval 
  # diagonal of the predicted cov is the variance of each prediction
  se.pred = 1.96 * sqrt(diag(cov.pred)) 
  poly.x = c(x.pred[1], x.pred, x.pred[length(x.pred)], rev(x.pred))
  poly.y = c((mu.pred - se.pred)[1], mu.pred + se.pred, 
             (mu.pred - se.pred)[length(x.pred)], rev(mu.pred - se.pred))
  polygon( poly.x, poly.y, col="#3182bd40", border="#3182bd50")
  
  # plot samples from the posterior
  # we have draws from an MvN using the learned mean and covariance of the function
  if (n.sample > 0) {
    for (i in 1:n.sample) {
      samp.y = mvrnorm(1, mu.pred, cov.pred)
      lines(x.pred, samp.y, col=col[i], lwd=1)
    }	
  }
  
  # plot the mean
  lines(x.pred, mu.pred, col="red", lwd=2)
  
  # plot the observations
  points(x, y, col="red", pch=19, cex=0.75)
}


# 1. fitting the GP with noise (y = k(x) + e)
# assuming noise is additive and iid
data("faithful", package = "datasets")
y = faithful$eruptions; x = faithful$waiting
summary(x)

sigma2e = 1
# x.test = matrix(seq(min(x), max(x), len = 200), nrow = 200, ncol = 1)
x.test = matrix(seq(0, 100, len = 200), nrow = 200, ncol = 1)
y.noisy = y + rnorm(length(x), mean = 0, sd = sqrt(sigma2e))

# solve gp and plot
gp = gp_solve(x, y.noisy, x.test, kfunc_exp, sigma2e )
gp_plot(x, y.noisy, x.test, gp[["mu"]], gp[["var"]], 
        main = "GP Regression with Sigma = 1", scale.se = F)


# 2. fitting the GP without noise (y=k(x))
gp = gp_solve(x, y, x.test, kfunc_exp)
gp_plot(x, y, x.test, gp[["mu"]], gp[["var"]], main="Gaussian Process Regression")
