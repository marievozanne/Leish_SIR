## SIMULATION STUDY -- Paper 1
## Model fit assessment -- SIR generated data, fit w/ MLR and SIR models

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
### Simulated data generated under SIR assumption -- first is wide form, second is long form, just weeks and categories
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/SIR_generated_simulation.RData")
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/lf_SIR_generated_simulation.RData")

## Load observed data and denote time and sample sizes
sim_sir <- sim_sir[sim_sir$Week!=0,]
y.sir.data <- sim_sir[sim_sir$Run==1,]
### Vector to designate weeks at which we have simulated data
week <- y.sir.data$Week
### Vector to designate simulated data sample sizes at each time point
n <- y.sir.data$mult_size
### Total population of Parnamirim
N <- 180000

## ----------------------------------------------------- ##
## ---------- Multinomial Logistic Regression ---------- ##
## Load MCMC object for SIR fit on multinomial generated data
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/sir_bmn.RData")

## Unique design matrix - all categories
W2 <- unique(model.matrix(Disease_Category ~ Week, long_sim_sir))

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

## Isolate MCMC estimates for multinomial logistic regression fit on SIR generated data
sir_mm_mcmc <- sir_bmn[,c(2,4,1,3)]
sir_mm_mcmc <- sir_mm_mcmc[1:nrow(sir_mm_mcmc) %% 10 == 0,] ## thin -- save every 10th iteration
sir_mm_mcmc <- sir_mm_mcmc[121:nrow(sir_mm_mcmc),] ## remove burn-in

## Sample from posterior
set.seed(123)
sir_mm_sample <- sample(1:nrow(sir_mm_mcmc), 5000, replace = TRUE)

## Make storage for predicted compartment proportions
sir_mm_P <- array(NA, dim=c(length(week), 3, length(sir_mm_sample)))

## Separate out sampled coefficient estimates that 
## correspond to S and I compartments
sir_mm_mcmc_coeff_S <- sir_mm_mcmc[sir_mm_sample,1:2]
sir_mm_mcmc_coeff_I <- sir_mm_mcmc[sir_mm_sample,3:4]

## Calculate predicted compartment proportions 
for (i in 1:length(sir_mm_sample)){
  sir.p <- P.fcn(W2, sir_mm_mcmc_coeff_S[i,], sir_mm_mcmc_coeff_I[i,])
  sir.p$PR <- 1-apply(sir.p, 1, sum)
  sir_mm_P[,,i] <- as.matrix(sir.p)
}

## Package up predictions for multinomial logistic fit
## using SIR generated data
sir_mm_P_pred_df <- data.frame(week, sir_mm_P[,,1])
names(sir_mm_P_pred_df) <- c("Week", "S", "I", "R")
for (i in 2:dim(sir_mm_P)[3]){
  temp3 <- data.frame(week,sir_mm_P[,,i])
  names(temp3) <- c("Week", "S", "I", "R")
  sir_mm_P_pred_df <- rbind(sir_mm_P_pred_df, temp3)
}
sir_mm_P_pred_df$Itr <- rep(1:5000, each=length(week))

## ------------------------------------------------- ##
## ------------ Compartmental SIR Model ------------ ##
source("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/derive_beta_parameters_from_quantiles.R")

geom.mean6weeks <- qgeom(c(0.05, 0.5, 0.95), 0.10)
geom.mean4weeks <- qgeom(c(0.05, 0.5, 0.95), 0.15)

strong_prior_param_calc <- beta.parms.from.quantiles(q=c(0.10, 0.15), p=c(0.05, 0.95), plot=F)

alpha_IR <- strong_prior_param_calc$a
beta_IR <- strong_prior_param_calc$b

## SIR fit on SIR generated data
# 1. Load S0, I0, gamma, beta, pi_IR from some saved file
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/SIR_simulation_fit123124_exp.RData")

# 2. Sample from posterior -- beta (pi_SI), S0, I0, pi_IR 
set.seed(123)
sir_cm_sample <- sample(335:dim(I_Array)[3], 5000, replace = TRUE)

wks <- 28

## Create storage for predicted proportions for SIR fit 
## on SIR generated data
sir_cm_P <- array(NA, dim=c(wks, 3, length(sir_cm_sample)))

for (idx in 1:length(sir_cm_sample)){
  S <- rep(NA, (wks+1))
  S[1] <- S0_Vector[sir_cm_sample[idx]]
  I <- rep(NA, (wks+1))
  I[1] <- I0_Vector[sir_cm_sample[idx]]
  
  I_SI <- rep(NA, wks)
  I_SI[1] <- 0
  R_IR <- rep(NA, wks)
  R_IR[1] <- 0
  
  for (i in 2:(wks+1)){
    I_SI[i] <- rbinom(1, S[i-1], pi_SI_Array[(i-1),,sir_cm_sample[idx]])
    R_IR[i] <- rbinom(1, I[i-1], pi_IR_Vector[sir_cm_sample[idx]])
    S[i] <- S[i-1]-I_SI[i]
    I[i] <- I[i-1]+I_SI[i]-R_IR[i]
  }
  
  I.aug <- I[-1]
  S.m <- S[-1]
  
  sir_cm_P[,,idx] <- cbind(S.m, I.aug, (N-S.m-I.aug))/N
}

sir_cm_P_array  <- array(data=NA, dim=c(wks, 4, length(sir_cm_sample)))
for (idx in 1:length(sir_cm_sample)){
  temp <- as.matrix(data.frame(1:wks,sir_cm_P[,,idx])) 
  sir_cm_P_array[,,idx] <- temp
}

## Package up predictions for SIR fit
## using SIR generated data 
sir_cm_P_pred_df <- sir_cm_P_array[,,1]
for (i in 2:length(sir_cm_sample)){
  sir_cm_P_pred_df <- rbind(sir_cm_P_pred_df, sir_cm_P_array[,,i])
}

sir_cm_P_pred_df <- data.frame(sir_cm_P_pred_df)
names(sir_cm_P_pred_df) <- c("Week", "S", "I", "R")
sir_cm_P_pred_df$Itr <- rep(1:5000, each=length(week))

## -------------------------------------------- ##
## ---------------- Make Plots ---------------- ##
### True values for generating multinomial data -- from 2005 survey
tp <- sum(0.510,0.019,0.350)
pi_S0 <- 0.510/tp ## re-weight
pi_I0 <- 0.019/tp
pi_R0 <- 0.350/tp

## True values for generating multinomial data with compartmental latent structure
sir.tpr <- y.sir.data[,1:4]
sir.tpr[,2:4] <- sir.tpr[,2:4]/N

## Plot predicted distributions and simulated data for 
## multinomial logistic fit and SIR generated data
(qq1 <- ggplot(sir_mm_P_pred_df, aes(x=as.factor(Week), y=S)) 
  + geom_violin(aes(x=as.factor(Week), y=S), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = S), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=S), data=sir.tpr, color="red", pch=8, size=3)
  # + xlab("Week") 
  + xlab("") 
  + ylab(expression(pi[S]))
  + scale_y_continuous(limits=c(0,0.72), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(qq2 <- ggplot(sir_mm_P_pred_df, aes(x=as.factor(Week), y=I)) 
  + geom_violin(aes(x=as.factor(Week), y=I), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = I), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=I), data=sir.tpr, color="red", pch=8, size=3)
  # + xlab("Week") 
  + xlab("") 
  + ylab(expression(pi[I]))
  + scale_y_continuous(limits=c(0,0.55), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(qq3 <- ggplot(sir_mm_P_pred_df, aes(x=as.factor(Week), y=R)) 
  + geom_violin(aes(x=as.factor(Week), y=R), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = R), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=R), data=sir.tpr, color="red", pch=8, size=3)
  + xlab("Week") + ylab(expression(pi[R]))
  + scale_y_continuous(limits=c(0,0.98), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
multiplot(qq1,qq2,qq3)

## Plot predicted distributions and simulated data for 
## SIR fit and SIR generated data
(q1 <- ggplot(sir_cm_P_pred_df, aes(x=as.factor(Week), y=S)) 
  + geom_violin(aes(x=as.factor(Week), y=S), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = S), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=S), data=sir.tpr, color="red", pch=8, size=3)
  # + xlab("Week") 
  + xlab("")
  + ylab(expression(pi[S]))
  + scale_y_continuous(limits=c(0,0.62), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(q2 <- ggplot(sir_cm_P_pred_df, aes(x=as.factor(Week), y=I)) 
  + geom_violin(aes(x=as.factor(Week), y=I), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = I), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=I), data=sir.tpr, color="red", pch=8, size=3)
  # + xlab("Week") 
  + xlab("")
  + ylab(expression(pi[I]))
  + scale_y_continuous(limits=c(0,0.55), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
(q3 <- ggplot(sir_cm_P_pred_df, aes(x=as.factor(Week), y=R)) 
  + geom_violin(aes(x=as.factor(Week), y=R), scale="width") 
  + geom_boxplot(aes(x = as.factor(Week), y = R), width = .1)
  + geom_point(mapping=aes(x=as.factor(Week), y=R), data=sir.tpr, color="red", pch=8, size=3)
  + xlab("Week") 
  + ylab(expression(pi[R]))
  + scale_y_continuous(limits=c(0,0.98), labels=scaleFUN)
  + theme_bw()
  + theme(axis.text=element_text(size=14), axis.title=element_text(size=18))
  + scale_x_discrete(breaks=seq(1, 28, 4)))
multiplot(q1, q2, q3)

multiplot(qq1, qq2, qq3, q1, q2, q3, cols=2)

## MAKE EPS SIMULATION PLOTS
setEPS()
postscript("M:/CPHS/Projects/Leish/Marie/Papers in progress/Multinomial LR vs Compartmental Model -- Leish/Journal of Applied Statistics/Ozanne_Figure_2a.eps")
multiplot(qq1,qq2,qq3,q1,q2,q3,cols = 2)
dev.off()

## MAKE SIR GENERATED SIMULATION PLOTS
# tiff("M:/CPHS/Projects/Leish/Papers in progress/Ozanne_Figure_2.tiff", 
#      height = 6, width = 10, units = 'in', type="windows", res=800,
#      compression="lzw")
# multiplot(qq1, qq2, qq3, q1, q2, q3, cols=2)
# dev.off()

## -------------------------------------------- ##
## ---------- Calculate Bayes Factor ---------- ## 

## Simulate from data posterior -- MLR
sir_mm_post <- data.frame(sim_data(100, sir_mm_P_pred_df))
names(sir_mm_post) <- c("Week", "S", "I", "R", "Itr")

## Simulate from data posterior -- SIR
sir_cm_post <- data.frame(sim_data(100, sir_cm_P_pred_df))
names(sir_cm_post) <- c("Week", "S", "I", "R", "Itr")

## Weights for using importance sampling -- MLR
w_sir_mm <- matrix(NA, ncol=ncol(sir_mm_mcmc), 
                   nrow=nrow(sir_mm_mcmc[sir_mm_sample,]))
for (i in 1:ncol(w_sir_mm)){
  w_sir_mm[,i] <- dnorm(sir_mm_mcmc[sir_mm_sample,i], 
                        mean=0, 
                        sd=sqrt(10))/(dnorm(sir_mm_mcmc[sir_mm_sample,i], 
                                            mean=mean(sir_mm_mcmc[sir_mm_sample,i]),
                                            sd=sd(sir_mm_mcmc[sir_mm_sample,i])))
}

## Normalize weights
W_sir_mm <- w_sir_mm[,1]*w_sir_mm[,2]*w_sir_mm[,3]*w_sir_mm[,4]
w_sir_mm <- W_sir_mm/sum(W_sir_mm)

## Weights for using importance sampling -- SIR
w_sir_cm <- array(NA, dim=c(length(sir_cm_sample), 5))

w_sir_cm[,1] <- dbinom(x=S0_Vector[sir_cm_sample], 
                       size=N, 
                       prob=pi_S0)/(dnorm(x=S0_Vector[sir_cm_sample], 
                                          mean=mean(S0_Vector[sir_cm_sample]), 
                                          sd=sd(S0_Vector[sir_cm_sample])))
w_sir_cm[,2] <- dbinom(x=I0_Vector[sir_cm_sample],
                       size=N-S0_Vector[sir_cm_sample], 
                       prob=pi_I0/(1-pi_S0))/(dnorm(x=I0_Vector[sir_cm_sample],
                                                    mean=mean(I0_Vector[sir_cm_sample]),
                                                    sd=sd(I0_Vector[sir_cm_sample])))
w_sir_cm[,3] <- dnorm(x=Beta_Matrix[1,sir_cm_sample], 
                      mean=0, sd=0.2)/(dnorm(x=Beta_Matrix[1,sir_cm_sample],
                                             mean=mean(Beta_Matrix[1,sir_cm_sample]),
                                             sd=sd(Beta_Matrix[1,sir_cm_sample])))
w_sir_cm[,4] <- dnorm(x=Beta_Matrix[2,sir_cm_sample], 
                      mean=0, sd=0.2)/(dnorm(x=Beta_Matrix[2,sir_cm_sample],
                                             mean=mean(Beta_Matrix[2,sir_cm_sample]),
                                             sd=sd(Beta_Matrix[2,sir_cm_sample])))
w_sir_cm[,5] <- dbeta(x=pi_IR_Vector[sir_cm_sample], 
                          shape1=alpha_IR, 
                          shape2=beta_IR)/(dnorm(x=pi_IR_Vector[sir_cm_sample],
                                                 mean=mean(pi_IR_Vector[sir_cm_sample]),
                                                 sd=sd(pi_IR_Vector[sir_cm_sample])))
# no gammas in simulation  

### Normalize weights
W_sir_cm <- w_sir_cm[,1]*w_sir_cm[,2]*w_sir_cm[,3]*w_sir_cm[,4]*w_sir_cm[,5]
w_sir_cm <- W_sir_cm/sum(W_sir_cm)

##### --- Calculate BFs using Euclidean norm --- #####
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
mean(c(AR_Beta, AR_I0, AR_S0, 1)) ## want to set eps such that prop_match is around this
## 0.4817

## Get as close as possible to smaller eps cutoff while ensuring at least 0.001 of indicators non-zero
mm_ind_sir <- eq.fcn(sir_mm_post, y.sir.data[,c(1,8:10)], eps=2600) 
mm_ind_sir$prop_match ##0.0012
cm_ind_sir <- eq.fcn(sir_cm_post, y.sir.data[,c(1,8:10)], eps=2250)
cm_ind_sir$prop_match ##0.4822

BF_sir <- matrix(mm_ind_sir$match, nrow=1)%*%matrix(w_sir_mm, ncol=1)/matrix(cm_ind_sir$match, nrow=1)%*%matrix(w_sir_cm, ncol=1)
1/BF_sir ##1448.77 -- indisputable evidence in favor of compartmental model