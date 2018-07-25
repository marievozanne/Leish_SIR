## Set working directory
setwd("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/")

## Load packaged and data
library(MCMCpack)
library(nnet)
load("human_new_2.RData")

## Frequentist Multivariate Analysis -- verify Lima paper results
human_new <- within(human_new, Area <- relevel(Area, ref = 3))
human_new <- within(human_new, Disease_category <- relevel(Disease_category, ref="A"))

fmn1 <- multinom(Lima_category ~ as.factor(Area) + Sex + Number_of_Residents + Family_Group_Income
                 + Time_Residing_In_Neighborhood_2 + as.factor(Cat_Age), data = human_new)
summary(fmn1)

## Frequentist Multivariate Analysis -- three categories
fmn2 <- multinom(Disease_category ~ Area + Sex + Number_of_Residents + Family_Group_Income
                 + Time_Residing_In_Neighborhood_2 + Cat_Age, data = human_new)
summary(fmn2)

fmn3 <- multinom(Disease_category ~ Area + Sex + Cat_Age + Week, data = human_new)
summary(fmn3)

## Bayesian Multivariate Analysis
bmn <- MCMCmnl(Disease_category ~ Area + Sex + Cat_Age + Week,
               mcmc.method="IndMH", B0=0, verbose=0,
               mcmc=100000, baseline="A", data = human_new)

summary(bmn)
