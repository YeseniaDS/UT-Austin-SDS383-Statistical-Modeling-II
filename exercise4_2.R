setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
library(ggplot2)
library(MASS)

# Waiting time between eruptions and the duration of the eruption 
# for the Old Faithful geyser in Yellowstone National Park, Wyoming, USA.
data("faithful", package = "datasets")
summary(faithful)

# plot the data points
ggplot(data = faithful) +
  geom_point(aes(waiting, eruptions), col = "blue")

# fit a linear model using lm
fit = lm(eruptions ~ waiting, data = faithful)
summary(fit)

# plot the fitted model
ggplot() +
  geom_point(data = faithful, aes(waiting, eruptions), col = "blue") +
  stat_smooth(data = faithful, aes(waiting, eruptions), method = "lm", se = FALSE)

# Guassian process (project into high dimension)
y = faithful$eruptions; x = faithful$waiting

# define the grid of values at which we will define the funtion
x.star = seq(min(x), max(x), len = 50)

# define phi
phi = cbind(1, x, x^2, x^3)

# define phi*
phi.star = cbind(1, x.star, x.star^2, x.star^3)

# set the inverse of posterior covariance as identity matrix
Sigma = diag(rep(1,4))
sigma2 = 1

# compute the posterior covariance 
A = (1/sigma2)*crossprod(phi) + solve(Sigma)
Ainv = solve(A)

# compute the estimated mean and covariance of the precitive distribution
m.star = (1/sigma2)* phi.star %*% Ainv %*% t(phi) %*% y
cov.star = phi.star %*% Ainv %*% t(phi.star)

# compute credible intervals 
interval = 1.96 * sqrt(diag(cov.star))
l1 = m.star - interval; l2 = m.star + interval

# plot the results
ggplot() +
  geom_point(data = faithful, aes(waiting, eruptions), col = "blue") +
  geom_line(aes(x.star, m.star)) +
  geom_ribbon(aes(x.star, ymin = l1, ymax = l2), alpha = 0.3)+
  ggtitle("Value Sequences Between [min(x), max(x)]") +
  theme(plot.title = element_text(hjust = 0.5))

# Or, we can sample from multivariate normal
n.samples = 50
values = matrix(rep(0, length(x.star) * n.samples), ncol = n.samples)

for (i in 1:n.samples) {
  values[,i] <- mvrnorm(1, m.star, cov.star)}

values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id ="x")

# plot the results
ggplot() +
  geom_point(data = faithful, aes(waiting, eruptions), col = "blue") +
  geom_line(data = values, aes(x = x, y = value, group = variable), col="grey") +
  geom_line(aes(x.star, m.star)) +
  ggtitle("Values Generated from MvN Distribution") +
  theme(plot.title = element_text(hjust = 0.5))
