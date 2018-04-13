setwd("D:/2018 UT Austin/R Code/Statistical Modeling II/section04")
library(readr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#---------------------------------------------------------------------
weather <- read.csv('weather.csv')
longitude <- weather$lon; latitude <- weather$lat
x <- seq(min(longitude) - 1 , max(longitude) + 1, 1)
y <- seq(min(latitude) - 1, max(latitude) + 1, 1)
summary(x); summary(y)
n1 = nrow(weather)
n2 = length(x) * length(y)
gp_data <- list(N1 = n1, D = 2, N2 = n2, y = weather$temperature,
                X1 = cbind(longitude, latitude), X2 = expand.grid(x, y))

#---------------------------------------------------------------------
# reduce the number of iterations to 500
fit <- stan("gp_fit_multi_inputs.stan", data = gp_data, iter = 500, chains = 4)
params <- extract(fit)
z <- matrix(colMeans(params$y2), length(x), length(y), byrow = F)

#---------------------------------------------------------------------
# break up to 20 intervals
intervals <- as.numeric(cut(weather$temperature, breaks = 20))

# plot mean function and contour
par(mar = c(2.1, 2.0, 2.1, 2.1))
pointcol <- colorRampPalette(c('red', 'blue'))(20)[intervals]
filled.contour(x, y, z, color.palette = colorRampPalette(c('yellow', 'darkgreen')),
               plot.axes = {points(longitude, latitude, pch = 17, col = pointcol)})

