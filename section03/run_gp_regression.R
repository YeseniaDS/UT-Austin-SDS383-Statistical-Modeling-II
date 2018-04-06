rm(list = ls())
library(readr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data(faithful)

gp_data<-read_rdump('D:/2018 UT Austin/R Code/Statistical Modeling II/faithful_data.R')

filename = "D:/2018 UT Austin/R Code/Statistical Modeling II/gp_regression.stan"

stan_code <-readChar(filename, file.info(filename)$size)

fit <- stan(model_code=stan_code,data=gp_data)

params<-extract(fit)

par(mfrow=c(1, 3))
alpha_breaks=10 * (0:50) / 50 - 5
hist(log(params$alpha), main="", xlab="log(alpha)", col="purple", yaxt='n')

beta_breaks=10 * (0:50) / 50 - 5
hist(log(params$rho), main="", xlab="log(rho)", col="purple", yaxt='n')

sigma_breaks=5 * (0:50) / 50
hist(params$sigma, main="", xlab="sigma", col="purple", yaxt='n')

readline(prompt="Press [enter] to continue")

probs = c(0.1,0.5,0.9)
cred <- sapply(1:length(gp_data$x_predict), function(n) quantile(params$y_predict[,n], probs=probs))

par(mfrow=c(1, 1))
plot(1, type="n",
     xlab="x", ylab="y", xlim=c(0, 120), ylim=c(-10, 10))
polygon(c(gp_data$x_predict, rev(gp_data$x_predict)), c(cred[1,], rev(cred[3,])),
        col = "grey", border = NA)
lines(gp_data$x_predict, cred[2,], lwd=2)

points(gp_data$x,gp_data$y)

