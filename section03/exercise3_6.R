rm(list = ls())
setwd("D:/2018 UT Austin/R Code/Statistical Modeling II")
set.seed(123)

# Data: numberof out of school suspensions
data <- read_csv('tea_discipline_oss.csv')
summary(data)

tea <- data %>% filter(ACTIONS >= 0) # Remove dropouts and NAs
summary(tea)
hist(tea$ACTIONS, breaks = 20, cex.main = 1)