## SIMULATION STUDY -- Paper 1
## Model fit assessment -- multinomial generated data, fit w/ MLR and SIR models

## Load libraries
library(ggplot2)
library(mc2d)
library(plyr)

## Load other functions
### Load SIR MCMC functions to calculate compartment proportions
source("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/SIR_MCMC_Functions.R")
### Load multiplot function for ggplot2
source("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/ggplot2_multiplot_fcn.R")
## Define scale function for plotting in ggplot2
scaleFUN <- function(x) sprintf("%.2f", x)

## Load simulated data
### Simulated data generated under multinomial assumption -- first is wide form, second is long form
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/Multinomial_generated_simulation.RData")
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/lf_Multinomial_generated_simulation.RData")

## Load observed data and denote time and sample sizes
y.mlr.data <- sim_multinomial[sim_multinomial$Run==1,]
### Vector to designate weeks at which we have simulated data
week <- y.mlr.data$Week
### Vector to designate simulated data sample sizes at each time point
n <- y.mlr.data$mult_size
### Total population of Parnamirim
N <- 180000

## ----------------------------------------------------- ##
## ---------- Multinomial Logistic Regression ---------- ##
## Load MCMC object for MLR fit on multinomial generated data
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/multinom_bmn.RData")

## Unique design matrix - all categories
W1 <- unique(model.matrix(Disease_Category ~ Week, long_sim_multinomial))

## Function for calculating P(S), P(I) for multinomial logistic regression
P.fcn <- function(X, b, g){
  eXb <- exp(X%*%b)
  eZg <- exp(X%*%g)
  PS <- eXb/(1+eXb)*(1/(1+eZg))/(1-(eXb/(1+eXb))*(eZg/(1+eZg)))
  PI <- (eZg/(1+eZg))*(1/(1+eXb))/(1-eXb/(1+eXb)*eZg/(1+eZg))
  return(data.frame(PS, PI))
}

## Function for generating from the posterior for the data
sim_data <- function(nn, P){
  set.seed(123)
  mult <- matrix(NA, nrow=nrow(P), ncol=ncol(P))
  for (i in 1:nrow(P)){
    mult[i,] <- c(P[i,1],rmultinom(n=1, size=nn, prob=P[i,2:4]),P[i,5])
  }
  mult
}

## Isolate MCMC estimates for multinomial logistic regression fit on multinomial generated data
mlr_bmn_mcmc <- multinom_bmn[,c(2,4,1,3)]
mlr_bmn_mcmc <- mlr_bmn_mcmc[1:nrow(mlr_bmn_mcmc) %% 10 == 0,] ## thin -- save every 30th iteration
mlr_bmn_mcmc <- mlr_bmn_mcmc[121:nrow(mlr_bmn_mcmc),] ## remove burn-in

## Sample from posterior
set.seed(123)
mlr_sample <- sample(1:nrow(mlr_bmn_mcmc), 5000, replace = TRUE)

## Make storage for predicted compartment proportions
mlr_P <- array(NA, dim=c(length(week), 3, length(mlr_sample)))

## Separate out sampled coefficient estimates that 
## correspond to S and I compartments
mlr_mcmc_coeff_S <- mlr_bmn_mcmc[mlr_sample,1:2]
mlr_mcmc_coeff_I <- mlr_bmn_mcmc[mlr_sample,3:4]

## Calculate predicted compartment proportions 
for (i in 1:length(mlr_sample)){
  mlr.p <- P.fcn(W1, mlr_mcmc_coeff_S[i,], mlr_mcmc_coeff_I[i,])
  mlr.p$PR <- 1-apply(mlr.p, 1, sum)
  mlr_P[,,i] <- as.matrix(mlr.p)
}

## Package up predictions for multinomial logistic fit
## using multinomial generated data
mlr_mm_P_pred_df <- data.frame(week, mlr_P[,,1])
names(mlr_mm_P_pred_df) <- c("Week", "S", "I", "R")
for (i in 2:dim(mlr_P)[3]){
  temp3 <- data.frame(week,mlr_P[,,i])
  names(temp3) <- c("Week", "S", "I", "R")
  mlr_mm_P_pred_df <- rbind(mlr_mm_P_pred_df, temp3)
}
mlr_mm_P_pred_df$Itr <- rep(1:5000, each=length(week))

## ------------------------------------------------- ##
## ------------ Compartmental SIR Model ------------ ##
source("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/derive_beta_parameters_from_quantiles.R")

geom.mean6weeks <- qgeom(c(0.05, 0.5, 0.95), 0.10)
geom.mean4weeks <- qgeom(c(0.05, 0.5, 0.95), 0.15)

strong_prior_param_calc <- beta.parms.from.quantiles(q=c(0.10, 0.15), p=c(0.05, 0.95), plot=F)

alpha_IR <- strong_prior_param_calc$a
beta_IR <- strong_prior_param_calc$b

## SIR fit on multinomial generated data
# 1. Load S0, I0, gamma, beta, pi_IR from some saved file
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/MLR_simulation_fit123124_exp.RData")

# 2. Sample from posterior -- beta (pi_SI), S0, I0, pi_IR 
set.seed(123)
sir_sample <- sample(335:dim(I_Array)[3], 5000, replace = TRUE)

wks <- 28

## Create storage for predicted proportions for SIR fit 
## on multinomial generated data
P.pred.mlr <- array(NA, dim=c(wks, 3, length(sir_sample)))

for (idx in 1:length(sir_sample)){
  S <- rep(NA, (wks+1))
  S[1] <- S0_Vector[sir_sample[idx]]
  I <- rep(NA, (wks+1))
  I[1] <- I0_Vector[sir_sample[idx]]
  
  I_SI <- rep(NA, wks)
  I_SI[1] <- 0
  R_IR <- rep(NA, wks)
  R_IR[1] <- 0
  
  for (i in 2:(wks+1)){
    I_SI[i] <- rbinom(1, S[i-1], pi_SI_Array[(i-1),,sir_sample[idx]])
    R_IR[i] <- rbinom(1, I[i-1], pi_IR_Vector[sir_sample[idx]])
    S[i] <- S[i-1]-I_SI[i]
    I[i] <- I[i-1]+I_SI[i]-R_IR[i]
  }
  
  I.aug <- I[-1]
  S.m <- S[-1]
  
  P.pred.mlr[,,idx] <- cbind(S.m, I.aug, (N-S.m-I.aug))/N
}

P.pred.mlr.array <- array(data=NA, dim=c(wks, 4, length(sir_sample)))
for (idx in 1:length(sir_sample)){
  temp <- as.matrix(data.frame(1:wks,P.pred.mlr[,,idx])) 
  P.pred.mlr.array[,,idx] <- temp
}

## Package up predictions for SIR fit
## using multinomial generated data
mlr_cm_P_pred_df <- P.pred.mlr.array[,,1]
for (i in 2:length(sir_sample)){
  mlr_cm_P_pred_df <- rbind(mlr_cm_P_pred_df, P.pred.mlr.array[,,i])
}

mlr_cm_P_pred_df <- data.frame(mlr_cm_P_pred_df)
names(mlr_cm_P_pred_df) <- c("Week", "S", "I", "R")
mlr_cm_P_pred_df$Itr <- rep(1:5000, each=length(week))

## -------------------------------------------- ##
## ---------------- Make Plots ---------------- ##
### True values for generating multinomial data -- from 2005 survey
tp <- sum(0.510,0.019,0.350)
pi_S0 <- 0.510/tp ## re-weight
pi_I0 <- 0.019/tp
pi_R0 <- 0.350/tp

## Make data frame with true constant proportions for each category and change names
mlr.tpr <- data.frame(c(1:28), c(rep(pi_S0, 28)), c(rep(pi_I0, 28)), c(rep(pi_R0, 28)))
names(mlr.tpr) <- c("Week", "S", "I", "R")

## Plot predicted distributions and simulated data for 
## multinomial logistic fit and multinomial generated data
(pp1 <- ggplot(mlr_mm_P_pred_df, aes(x=as.factor(Week), y=S)) 
  + geom_violin(aes(x=as.factor(Week), y=S), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = S), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=S), data=mlr.tpr, color="red", pch=8, size=3)
  # + xlab("Week")
  + xlab("")
  + ylab(expression(pi[S]))
  + scale_y_continuous(limits=c(0,0.72), labels = scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(pp2 <- ggplot(mlr_mm_P_pred_df, aes(x=as.factor(Week), y=I)) 
  + geom_violin(aes(x=as.factor(Week), y=I), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = I), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=I), data=mlr.tpr, color="red", pch=8, size=3)
  # + xlab("Week")
  + xlab("") 
  + ylab(expression(pi[I]))
  + scale_y_continuous(limits=c(0,0.1), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(pp3 <- ggplot(mlr_mm_P_pred_df, aes(x=as.factor(Week), y=R)) 
  + geom_violin(aes(x=as.factor(Week), y=R), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = R), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=R), data=mlr.tpr, color="red", pch=8, size=3)
  + xlab("Week") + ylab(expression(pi[R]))
  + scale_y_continuous(limits=c(0,0.75), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
multiplot(pp1, pp2, pp3)

## Plot predicted distributions and simulated data for 
## SIR fit and multinomial generated data
(p1 <- ggplot(mlr_cm_P_pred_df, aes(x=as.factor(Week), y=S)) 
  + geom_violin(aes(x=as.factor(Week), y=S), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = S), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=S), data=mlr.tpr, color="red", pch=8, size=3)
  + xlab("Week") + ylab(expression(pi[S]))
  + scale_y_continuous(limits=c(0,0.65), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(p2 <- ggplot(mlr_cm_P_pred_df, aes(x=as.factor(Week), y=I)) 
  + geom_violin(aes(x=as.factor(Week), y=I), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = I), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=I), data=mlr.tpr, color="red", pch=8, size=3)
  + xlab("Week") + ylab(expression(pi[I]))
  + scale_y_continuous(limits=c(0,0.1), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(p3 <- ggplot(mlr_cm_P_pred_df, aes(x=as.factor(Week), y=R)) 
  + geom_violin(aes(x=as.factor(Week), y=R), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = R), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=R), data=mlr.tpr, color="red", pch=8, size=3)
  + xlab("Week") + ylab(expression(pi[R]))
  + scale_y_continuous(limits=c(0,0.75), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
multiplot(p1,p2,p3)

multiplot(pp1,pp2,pp3,p1,p2,p3,cols = 2)

## MAKE EPS SIMULATION PLOTS
setEPS()
postscript("M:/CPHS/Projects/Leish/Marie/Papers in progress/Multinomial LR vs Compartmental Model -- Leish/Journal of Applied Statistics/Ozanne_Figure_1a.eps")
multiplot(pp1,pp2,pp3,p1,p2,p3,cols = 2)
dev.off()

## MAKE MM GENERATED SIMULATION PLOTS
# tiff("M:/CPHS/Projects/Leish/Papers in progress/Ozanne_Figure_1.tiff", 
#      height = 6, width = 10, units = 'in', type="windows", res=800,
#      compression="lzw")
# multiplot(pp1, pp2, pp3, p1, p2, p3, cols=2)
# dev.off()
## -------------------------------------------- ##
## ---------- Calculate Bayes Factor ---------- ## 

## Simulate from data posterior -- MLR
bmn_mlr_post <- data.frame(sim_data(100, mlr_mm_P_pred_df))
names(bmn_mlr_post) <- c("Week", "S", "I", "R", "Itr")

## Simulate from data posterior -- SIR
cm_mlr_post <- data.frame(sim_data(100, mlr_cm_P_pred_df))
names(cm_mlr_post) <- c("Week", "S", "I", "R", "Itr")

## Weights for using importance sampling -- MLR
w_bmn_mlr <- matrix(NA, ncol=ncol(mlr_bmn_mcmc), 
                    nrow=nrow(mlr_bmn_mcmc[mlr_sample,]))
for (i in 1:ncol(w_bmn_mlr)){
  w_bmn_mlr[,i] <- dnorm(mlr_bmn_mcmc[mlr_sample,i], 
                         mean=0, 
                         sd=sqrt(10))/(dnorm(mlr_bmn_mcmc[mlr_sample,i], 
                                             mean=mean(mlr_bmn_mcmc[mlr_sample,i]),
                                             sd=sd(mlr_bmn_mcmc[mlr_sample,i])))
}
W_bmn_mlr <- w_bmn_mlr[,1]*w_bmn_mlr[,2]*w_bmn_mlr[,3]*w_bmn_mlr[,4]
w_bmn_mlr <- W_bmn_mlr/sum(W_bmn_mlr)


## Weights for using importance sampling -- SIR
w_cm_mlr <- array(NA, dim=c(length(mlr_sample), 5))

w_cm_mlr[,1] <- dbinom(x=S0_Vector[mlr_sample], 
                       size=N, 
                       prob=pi_S0)/(dnorm(x=S0_Vector[mlr_sample], 
                                          mean=mean(S0_Vector[mlr_sample]), 
                                          sd=sd(S0_Vector[mlr_sample])))
w_cm_mlr[,2] <- dbinom(x=I0_Vector[mlr_sample],
                       size=N-S0_Vector[mlr_sample], 
                       prob=pi_I0/(1-pi_S0))/(dnorm(x=I0_Vector[mlr_sample],
                                                    mean=mean(I0_Vector[mlr_sample]),
                                                    sd=sd(I0_Vector[mlr_sample])))
w_cm_mlr[,3] <- dnorm(x=Beta_Matrix[1,mlr_sample], 
                      mean=0, sd=0.2)/(dnorm(x=Beta_Matrix[1,mlr_sample],
                                             mean=mean(Beta_Matrix[1,mlr_sample]),
                                             sd=sd(Beta_Matrix[1,mlr_sample])))
w_cm_mlr[,4] <- dnorm(x=Beta_Matrix[2,mlr_sample], 
                      mean=0, sd=0.2)/(dnorm(x=Beta_Matrix[2,mlr_sample],
                                             mean=mean(Beta_Matrix[2,mlr_sample]),
                                             sd=sd(Beta_Matrix[2,mlr_sample])))
w_cm_mlr[,5] <- dbeta(x=pi_IR_Vector[mlr_sample], 
                      shape1=alpha_IR, 
                      shape2=beta_IR)/(dnorm(x=pi_IR_Vector[mlr_sample],
                                             mean=mean(pi_IR_Vector[mlr_sample]),
                                             sd=sd(pi_IR_Vector[mlr_sample])))
  # no gammas in simulation  

### Normalize weights
W_cm_mlr <- w_cm_mlr[,1]*w_cm_mlr[,2]*w_cm_mlr[,3]*w_cm_mlr[,4]*w_cm_mlr[,5]
w_cm_mlr <- W_cm_mlr/sum(W_cm_mlr)

##### --- Calculate BFs using spectral norm for distance --- #####
### Calculate indicators -- do resampled data and real data agree?
### Function to check equality of samples (Euclidean norm):
eq.fcn <- function(post_data, y_data, eps){
  match <- rep(0, length=max(post_data$Itr))
  for (i in unique(post_data$Itr)){
    if (sum(colSums((post_data[post_data$Itr==i,2:4]-y_data[,2:4])^2)) < eps){
      match[i] <- 1
    }
  }
  return(list(match=match, prop_match=mean(match)))
}


## Calculation include AR of 1 for Gibbs sampling component
mean(c(AR_Beta, AR_I0, AR_S0, 1))
## 0.6187
mean(c(AR_Beta, AR_I0, AR_S0))
##0.4917

bmn_ind_mlr <- eq.fcn(bmn_mlr_post, y.mlr.data, eps=2600) ##2430; 0.4916
bmn_ind_mlr$prop_match ##0.6106
cm_ind_mlr  <- eq.fcn(cm_mlr_post,  y.mlr.data, eps=2600)
cm_ind_mlr$prop_match  ##0.3036

BF_multinom <- matrix(bmn_ind_mlr$match, nrow=1)%*%matrix(w_bmn_mlr, ncol=1)/matrix(cm_ind_mlr$match, nrow=1)%*%matrix(w_cm_mlr, ncol=1)
BF_multinom   ## 1.4605 -- no compelling reason to prefer one over the other
1/BF_multinom ## 0.6847