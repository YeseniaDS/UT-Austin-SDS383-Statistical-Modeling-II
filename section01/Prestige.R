## Statistical Modeling II exercise 1.25 
# load the library
library(car)

# load the data
prestige = data(Prestige)

# get rid of categorical variables
prestige = Prestige[,1:4]
head(prestige)

# split data
y = prestige[,4] # using income as dependent variable
x = as.matrix(prestige[,1:3]) # using the three non-categorical predictors

# add an intercept
x = cbind(1,x)

# compute the estimator
betahat = solve(crossprod(x)) %*% crossprod(x,y)

# compute the estimated standard error
e = y - x %*% betahat # residuals
betacov = sqrt(crossprod(e)/(nrow(x)-ncol(x)))
print(betacov) # 7.8465

# now compare to lm
fit = lm(y ~ x-1) # the 'minus 1' notation says not to fit an intercept 

summary(fit)
betacovlm = vcov(fit)
sqrt(diag(betacovlm))

# fetch the built in standard deviations of the residuals
sd(fit$residuals) # 7.7291
     