# part of the code credits to Aleksandar Dimitriev
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section04")
library(readr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#---------------------------------------------------------------------
data(faithful)
gp_data <- list(N = 10, N_predict = 121, x = faithful$waiting[1:10],
                y = faithful$eruptions[1:10], x_predict = seq(0,120, 1))
fit <- stan("gp_regression.stan", data = gp_data)
params <- extract(fit)

probs <- c(0.025, 0.5, 0.975)
qa <- quantile(params$alpha, probs)
qr <- quantile(params$rho,   probs)
qs <- quantile(params$sigma, probs)

#---------------------------------------------------------------------
# plot the histogram of alpha, rho, and sigma
par(mar = c(4.1, 4.0, 2.1, 2.1))
par(mfrow = c(1, 3))
hist(params$alpha, yaxt = 'n', breaks = 100, cex.main = 1,
     main = sprintf("alpha 95%% CI: (%.2f, %.2f)", qa[1], qa[3]))
hist(params$rho, yaxt = 'n', breaks = 100, cex.main = 1,
     main = sprintf("rho 95%% CI: (%.2f, %.2f)", qr[1], qr[3]))
hist(params$sigma, yaxt = 'n', breaks = 100, cex.main = 1,
     main = sprintf("sigma 95%% CI: (%.2f, %.2f)",qs[1], qs[3]))

# traceplot
traceplot(fit, par = c("rho", "alpha", "sigma"))

#---------------------------------------------------------------------
cred <- sapply(1:length(gp_data$x_predict), 
               function(n) quantile(params$y_predict[,n], probs = probs))
# plot 95% credible interval
par(mfrow = c(1, 1))
plot(1, type = "n", xlab = "x", ylab = "y", xlim = c(0, 120), 
     ylim = c(-10, 10), main = "Stan Gaussian Process", cex.main = 1)
polygon(c(rev(gp_data$x_predict), gp_data$x_predict), c(rev(cred[3,]), cred[1,]),
        col = "#3182bd40", border = "#3182bd50")
lines(gp_data$x_predict, cred[2,], lwd = 2)
points(gp_data$x,gp_data$y, pch = 19, col = "red")


