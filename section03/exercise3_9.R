library(readr)
library(rstan)
tea_discipline_oss <- read.csv("/home/xinjie/Dropbox/courses/2018_Spring/STAT_MODEL/sds383d/code/ex3/tea_discipline_oss.csv") 
View(tea_discipline_oss)

uncensored_data = subset(tea_discipline_oss,ACTIONS>0)
tea <-data.frame(x=as.numeric(as.character(uncensored_data$GRADE)),y=uncensored_data$ACTIONS)
#tea <-data.frame(y=uncensored_data$ACTIONS)
tea$intercept =1
tea<-as.list(tea)
tea$N<-nrow(uncensored_data)
fileName <- "poisson.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
resStan<-stan(model_code=stan_code,data=tea,chains=3,iter=1500,warmup=500,thin=10)
traceplot(resStan, pars = c("beta"), inc_warmup = FALSE) #set inc_warmup = TRUE to see burn in
#                mean  se_mean       sd         2.5%          25%          50%          75%        97.5% n_eff     Rhat
#beta[1]   2.389668e+00 0.000357 0.004925 2.380165e+00 2.386489e+00 2.389721e+00 2.392995e+00 2.399572e+00   191 0.999809
#beta[2]   5.022300e-02 0.000044 0.000633 4.890300e-02 4.982300e-02 5.022300e-02 5.065700e-02 5.143500e-02   207 0.996889

# include more variables
tea$sex = as.numeric(uncensored_data$SEXX == 'MALE')
fileName <- "poisson.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
resStan<-stan(model_code=stan_code,data=tea,chains=3,iter=1500,warmup=500,thin=10)
traceplot(resStan, pars = c("beta"), inc_warmup = FALSE) #set inc_warmup = TRUE to see burn in

#                  mean  se_mean       sd         2.5%          25%          50%
#beta[1]   2.306000e+00 0.000363 0.005807 2.295024e+00 2.302043e+00 2.305959e+00
#beta[2]   5.193500e-02 0.000038 0.000619 5.073100e-02 5.151700e-02 5.188200e-02
#beta[3]   9.821100e-02 0.000228 0.003863 9.092200e-02 9.570600e-02 9.803000e-02
#               75%        97.5%       n_eff     Rhat
#beta[1]   2.310049e+00 2.317261e+00   255 1.004584
#beta[2]   5.237900e-02 5.312200e-02   270 1.005256
#beta[3]   1.008720e-01 1.057470e-01   288 0.992144