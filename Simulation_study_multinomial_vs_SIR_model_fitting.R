## Simulation study model fitting
## Multinomial-generated and SIR-generated simulation data 
## to be fit using bayesian multinomial logistic regression
## and compartmental model

## Load packages
library(MCMCpack)
library(reshape2)

## --- Multinomial logistic regression --- 

# with SIR compartmental model generated data ---
# Reformat data -- long form
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/SIR_generated_simulation.RData")
m_sim_sir <- melt(sim_sir[sim_sir$Run==1,c(1,8:10)], 
                  id.vars=c("Week"), variable.name=c("Disease_Category"),
                  value.name=c("Count"))
m_sim_sir <- m_sim_sir[m_sim_sir$Week != 0,]

long_sim_sir <- NULL
for (i in 1:nrow(m_sim_sir)){
  if (m_sim_sir$Count[i] > 0){
    temp <- data.frame(matrix(unlist(rep(m_sim_sir[i,1], m_sim_sir$Count[i])), ncol=1))
    names(temp) <- c("Week")
    temp$Disease_Category <- rep(m_sim_sir[i,2], m_sim_sir$Count[i])
    long_sim_sir <- rbind(long_sim_sir, temp)
  }
  else {
    long_sim_sir <- long_sim_sir
  }
}
save(long_sim_sir, file="M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/lf_SIR_generated_simulation.RData")

# Fit multinomial logistic regression model -- SIR generated data

ptm <- proc.time()
sir_bmn <- MCMCmnl(Disease_Category ~ Week, 
                   mcmc.method="IndMH", B0=0.1, verbose=0,
                   mcmc=150000, baseline="Y_R", data=long_sim_sir,
                   seed=123124)
proc.time() - ptm

summary(sir_bmn)
save(sir_bmn, file="M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/sir_bmn.RData")

# with multinomial model generated data ---
# Reformat data -- long form
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/Multinomial_generated_simulation.RData")
m_sim_multinomial <- melt(sim_multinomial[sim_multinomial$Run==1, c(1:4)],
                          id.vars=c("Week"), variable.name=c("Disease_Category"),
                          value.name=c("Count"))

long_sim_multinomial <- NULL
for (i in 1:nrow(m_sim_multinomial)){
  if (m_sim_multinomial$Count[i] > 0){
    temp <- data.frame(matrix(unlist(rep(m_sim_multinomial[i,1], m_sim_multinomial$Count[i])), ncol=1))
    names(temp) <- c("Week")
    temp$Disease_Category <- rep(m_sim_multinomial[i,2], m_sim_multinomial$Count[i])
    long_sim_multinomial <- rbind(long_sim_multinomial, temp)
  }
  else{
    long_sim_multinomial <- long_sim_multinomial
  }
}
save(long_sim_multinomial, file="M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/lf_Multinomial_generated_simulation.RData")

multinom_bmn <- MCMCmnl(Disease_Category ~ Week,
                        mcmc.method = "IndMH", B0=0.1, verbose=0,
                        mcmc=150000, baseline = "Y_R", data=long_sim_multinomial,
                        seed=123124)
summary(multinom_bmn)
save(multinom_bmn, file="M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/multinom_bmn.RData")

### ----- Compartmenal Model Fitting ----- ### 
runMCMCChain <- function(seed){
  set.seed(seed);
  setwd("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/")
  # library(coda)
  library(mc2d)
  library(plyr)
  
  ## ------- Deterministic Update Functions ------- 
  
  calculate_S <- function(S0, I_SI){
    s <- rbind(S0,
               (matrix(S0, 
                       nrow = nrow(I_SI), 
                       ncol = length(S0), byrow = TRUE) - 
                  (apply(I_SI, 2, cumsum))))[1:nrow(I_SI),]
    s <- as.matrix(s)
    s
  }
  
  calculate_I <- function(I0, I_SI, R_IR){
    i <- matrix(rbind(I0,
                      (matrix(I0,
                              nrow = nrow(R_IR),
                              ncol = length(I0), byrow = TRUE) +
                         (apply(I_SI, 2, cumsum) - 
                         (apply(R_IR, 2, cumsum)))))[1:nrow(R_IR),])
    i <- as.matrix(i)
    i
  }
  
  calculate_R <- function(R0, R_IR){
    r <- matrix(rbind(R0,
                      (matrix(R0,
                              nrow = nrow(R_IR),
                              ncol = length(R0), byrow = TRUE) + 
                         (apply(R_IR, 2, cumsum)))))[1:nrow(R_IR),]
    r <- as.matrix(r)
    r
  }
  
  ## ------- Set up of probability calculations dependent on time -------
  
  ## --- Sandfly components ---
  betaparms <- function(mu, var){
    a <- mu/(1-mu)*((1-var*(1-mu))/var)
    b <- (1-var*(1-mu))/var
    if (a > 0 && b > 0){
      return(c(a,b))
    }
    else{
      paste("That's silly!")
    }
  }
  
  ## --- Exponential transition probability ---
  calculate_pi_SI_exp <- function(X, Beta, N, I0, I_SI, R_IR, Itr){
    I <- calculate_I(I0, I_SI, R_IR)
    if (any(I < 0)){return(-Inf)}
    
    eta <- as.matrix(X[-Itr,])%*%matrix(Beta, nrow=ncol(X))
    
    pi_SI.fcn <- function(i, n, eta){
      p <- 1-exp(-i/n*exp(eta))
      p
    }
    matrix(mapply(pi_SI.fcn, I, N, eta),
           nrow = nrow(I),
           ncol = 1, byrow = TRUE)
  }
  
  ## --- Beta-Poisson transition probability ---
  calculate_pi_SI_BP <- function(X, Beta, Brate, I0, I_SI, R_IR, 
                                 N, delta, pinf_dogs, Itr){
    I <- calculate_I(I0, I_SI, R_IR)
    if (any(I < 0)){return(-Inf)}
    
    eta <- as.matrix(X[-Itr,])%*%matrix(Beta, nrow=ncol(X))
    mu_Binf <- (delta*I/N+(1-delta)*pinf_dogs)*exp(eta)
    var_Binf <- 1 ## Think about this more.
    inf_parms <- t(mapply(betaparms, mu_Binf, var_Binf))
    
    pi_SI.fcn <- function(rate, alpha, beta){
      p <- 1-(1+rate/beta)^(-alpha)
      p
    }
    matrix(mapply(pi_SI.fcn, Brate, inf_parms[,1], inf_parms[,2]),
           nrow=nrow(I),
           ncol=1, byrow=TRUE)
  }
  
  ## ------- Transition probability prior distributions -------
  
  ## Specify strong priors for shape and scale parameters for transition probabilities (beta prior)
  Strong_prior_fcn <- function(mean, ESS){
    alpha <- ESS*mean
    beta <- ESS*(1-mean)
    c(alpha, beta)
  }
  
  ## ------- Full conditional functions -------
  
  ### --- Transition probabilities ---
  ## Beta_FC - choice of exponential or beta-Poisson for pi_SI
  calculate_pi_SI <- function(X, Beta, I0, I_SI, R_IR, N, Brate, delta, pinf_dogs, Itr, 
                              fclform = c("exponential", "betaPoisson")){
    fclform <- match.arg(fclform)
    switch(fclform,
           exponential=calculate_pi_SI_exp(X, Beta, N, I0, I_SI, R_IR, Itr),
           betaPoisson=calculate_pi_SI_BP(X, Beta, Brate, I0, I_SI, R_IR,
                                          N, delta, pinf_dogs, Itr))
  }
  
  ## Beta_FC
  Beta_FC <- function(Beta, beta_mean, beta_sd, N, X, Y, 
                      S0, I0, I_SI, R_IR, pi_S0, pi_I0,
                      Brate, delta, pinf_dogs, Itr, time, fclform){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    if (S0 + I0 > N){return(-Inf)}
    
    proposed_pi_SI <- calculate_pi_SI(X, Beta, I0, I_SI, R_IR,
                                      N, Brate, delta, pinf_dogs, Itr, fclform)
    
    if (any(proposed_pi_SI > 1)){ 
      return(-Inf)
    }

    # Likelihood component:
    lik.data <- sum(dmultinomial(Y, prob = cbind(S[time,]/N, I[time,]/N, 
                                                 1-(S[time,]+I[time,])/N), log=TRUE))
    lik.ic <- (dbinom(S0, N, pi_S0, log=TRUE) + 
               dbinom(I0, N-S0, pi_I0/sum(c(pi_I0,1-pi_S0-pi_I0)), log=TRUE))
    lik.tc <- sum(dbinom(I_SI, S, proposed_pi_SI, log = TRUE))
    prior <- sum(dnorm(Beta, mean = c(0,-1), sd = c(1,1), log = TRUE))
    return(lik.data+lik.ic+lik.tc+prior)
  }
  
  ## pi_IR_FC
  pi_IR_FC <- function(pi_IR, S0, I0, I_SI, R_IR, pi_I0,
                       N, Y, alpha_IR, beta_IR, time){
    if (pi_IR < 0){
      return(-Inf)
    }
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    if (S0 + I0 > N){return(-Inf)}
    
    # Likelihood Component:
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I[time,]/N, 1-(S[time,]+I[time,])/N), log = TRUE))
    lik.ic <- dbinom(I0, N-S0, pi_I0/sum(c(pi_I0,1-pi_S0-pi_I0)), log=TRUE)
    lik.tc <- sum(dbinom(R_IR, I, pi_IR, log = TRUE))
    prior <- sum(dbeta(pi_IR, alpha_IR, beta_IR, log = TRUE))
    return(lik.data+lik.ic+lik.tc+prior)
  }
  
  ### --- Transition compartments --- 
  ## I_SI_FC
  I_SI_FC <- function(I_SI, R_IR, N, X, Y, Beta, S0, I0, pi_IR, pi_S0, pi_I0,
                      Brate, delta, pinf_dogs, time, Itr, fclform){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    if (S0 + I0 > N){return(-Inf)}
    
    pi_SI <- calculate_pi_SI(X, Beta, I0, I_SI, R_IR, N, Brate, delta, 
                             pinf_dogs, Itr, fclform)
    
    if (any(pi_SI > 1)){
      return(-Inf)
    }
    
    #Likelihood component
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I[time,]/N, 1-(S[time,]+I[time,])/N), log = TRUE))
    lik.ic <- (dbinom(S0, N, pi_S0, log=TRUE) + 
                 dbinom(I0, N-S0, pi_I0/sum(c(pi_I0,1-pi_S0-pi_I0)), log=TRUE))
    lik.tc <- sum(dbinom(I_SI, S, pi_SI, log = TRUE)) +
      sum(dbinom(R_IR, I, pi_IR, log = TRUE))
    return(lik.data+lik.ic+lik.tc)
  }
  
  ## R_IR_FC 
  R_IR_FC <- function(R_IR, I_SI, S0, I0, Y, N, pi_IR, pi_S0, time){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    if (S0 + I0 > N){return(-Inf)}
    
    #Likelihood component
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I[time,]/N, 
                                               1-(S[time,]+I[time,])/N), log = TRUE))
    lik.ic <- dbinom(I0, N-S0, pi_I0/sum(c(pi_I0,1-pi_S0-pi_I0)), log=TRUE)
    lik.tc <- sum(dbinom(R_IR, I, pi_IR, log = TRUE))
    return(lik.data+lik.ic+lik.tc)
  }
  
  ### --- Initial compartments ---
  S0_FC <- function(S0, N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X, Brate, 
                    delta, pinf_dogs, Itr, time, fclform){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    
    pi_SI <- calculate_pi_SI(X, Beta, I0, I_SI, R_IR, N, Brate, delta, pinf_dogs, Itr, fclform)
    if (any(pi_SI > 1)){return(-Inf)}
    
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I[time,]/N, 
                                               1-(S[time,]+I[time,])/N), log = TRUE))
    lik.ic <- sum(dbinom(S0, N, pi_S0, log=TRUE))
    lik.tc <- sum(dbinom(I_SI, S, pi_SI, log=TRUE))
    return(lik.data+lik.ic+lik.tc)
  }
  
  I0_FC <- function(I0, S0, N, I_SI, R_IR, pi_I0, pi_IR, time){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I[time,]/N, 
                                               1-(S[time,]+I[time,])/N), log = TRUE))
    lik.ic <- sum(dbinom(I0, N, pi_I0, log=TRUE))
    lik.tc <- sum(dbinom(R_IR, I, pi_IR, log=TRUE))
    return(lik.data+lik.ic+lik.tc)
  }
  
  ### ------- Proposal Functions -------
  
  Propose_Coef <- function(Coef, proposal_Mean, proposal_SD){
    perturbation <- rnorm(length(Coef), proposal_Mean, proposal_SD)
    proposed_Coef <- Coef + perturbation
    proposed_Coef  
  }
  
  Propose_Initial_Compartment <- function(N, Initial_Compartment, initial_prob){
    proposal <- rbinom(1, N, initial_prob)
    
    prob_proposed <- dbinom(proposal, N, initial_prob, log=TRUE)
    prob_current <- dbinom(Initial_Compartment, N, initial_prob, log=TRUE)
    
    list(Initial_Compartment=proposal,
         prob_current=prob_current,
         prob_proposed=prob_proposed)
  }
  
  Propose_New_Transition_Probability <- function(Transition_Compartment,
                                                 From_Compartment, 
                                                 alpha, beta){
    proposed_Transition_Probability <- rbeta(1, alpha+sum(Transition_Compartment), 
                                             beta+sum(From_Compartment-Transition_Compartment))
    proposed_Transition_Probability
  }
  
  Propose_New_Transition_Better_I_Promise_Sort_Of <- function(
    Transition_Compartment, 
    From_Compartment,
    proposal_proportion){
    idx <- sample(1:nrow(Transition_Compartment), 
                  size = floor(proposal_proportion*nrow(Transition_Compartment)), 
                  replace = FALSE)
    
    m <- ncol(Transition_Compartment)
    
    proposal <- Transition_Compartment
    
    proposalProb <- ifelse((From_Compartment[idx,] == 0), 
                           0, 
                           Transition_Compartment[idx,]/From_Compartment[idx,])
    if (any(Transition_Compartment > From_Compartment)){
      Transition_Compartment <<- Transition_Compartment
      From_Compartment <<- From_Compartment
      proposalProb <<- proposalProb
      print(Transition_Compartment)
      print("From:")
      print(From_Compartment)
      print(which(Transition_Compartment > From_Compartment))
      print(proposalProb)
      print(which(proposalProb > 1))
      stop("Invalid transition state: transition compartment is too greedy")
    }
    
    if (any(Transition_Compartment<0)){
      Transition_Compartment <<- Transition_Compartment
      From_Compartment <<- From_Compartment
      proposalProb <<- proposalProb
      print(Transition_Compartment)
      print("From:")
      print(From_Compartment)
      print(which(Transition_Compartment > From_Compartment))
      print(proposalProb)
      print(which(proposalProb > 1))
      stop("Invalid transition state: proposal probability is negative")
    }
    proposal[idx,] <- rbinom(rep(1, length(idx)),
                             size = From_Compartment[idx,],
                             prob = proposalProb)
    prob_proposed <- sum(dbinom(proposal[idx,],
                                size = From_Compartment[idx,],
                                prob = proposalProb, log = TRUE))
    
    #TransitionCompartment - proposal
    proposedFrom <- From_Compartment[idx,] #From_Compartment[idx,] + (Transition_Compartment[idx,] - proposal[idx,])
    reverseProposalProb <- ifelse(proposedFrom == 0, 
                                  0, 
                                  proposal[idx,]/proposedFrom)
    prob_current <-  sum(dbinom(Transition_Compartment[idx,],
                                size = proposedFrom,
                                prob = reverseProposalProb, log = TRUE))
    
    list(compartment=proposal, 
         prob_current = prob_current,
         prob_proposed = prob_proposed)
  }
  
  
  ### ------- Draw Functions -------
  
  # Human model
  Draw_Beta <- function(Beta, beta_mean, beta_sd, N, X, Y, S0, I0, I_SI, R_IR, pi_S0, pi_I0,
                        Brate, delta, pinf_dogs, proposal_Mean, proposal_SD, Itr, time, fclform){
    Beta_proposed <- Propose_Coef(Beta, proposal_Mean, proposal_SD)
    fc_0 <- Beta_FC(Beta,          beta_mean, beta_sd, N, X, Y, S0, I0, I_SI, R_IR, pi_S0, pi_I0, 
                    Brate, delta, pinf_dogs, Itr, time, fclform)
    fc_1 <- Beta_FC(Beta_proposed, beta_mean, beta_sd, N, X, Y, S0, I0, I_SI, R_IR, pi_S0, pi_I0,
                    Brate, delta, pinf_dogs, Itr, time, fclform)
    Rat <- fc_1 - fc_0
    
    if (log(runif(1)) < min(Rat,0) && !is.nan(Rat)){
      Beta <- Beta_proposed
      AR_Beta <- 1
    }
    else{
      Beta <- Beta
      AR_Beta <- 0
    }
    return(list(Beta, AR_Beta))
  }
  
  checkNA <- function(obj, name)
  {
    if (is.na(obj)){
      print(paste(name,"was NA:")) 
      print(obj)
    }
  }
  
  Draw_pi_IR <- function(pi_IR, N, Y, S0, I0, I_SI, R_IR, pi_I0, alpha_IR, beta_IR, time){
    I <- calculate_I(I0, I_SI, R_IR)
    
    pi_IR_proposed <- Propose_New_Transition_Probability(R_IR, I, alpha_IR, beta_IR)
    pi_IR_proposed
  }
  
  Draw_I_SI <- function(I_SI, R_IR, N, X, Y, Beta, S0, I0, pi_IR, pi_S0, pi_I0,
                        Brate, delta, pinf_dogs, time, pro_prop, fclform){
    S <- calculate_S(S0, I_SI)
    # if (any(S < 0)){return(list(-Inf, 0))}
    proposed_Trans <- Propose_New_Transition_Better_I_Promise_Sort_Of(
      I_SI, S, proposal_proportion = pro_prop)
    if (any(proposed_Trans$Compartment > S)){
      print("Propose_New_Transition_Better_I_Promise_Sort_Of is not doing its job.")
    }
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    I_SI_proposed <- proposed_Trans$compartment
    
    fc_0 <- I_SI_FC(I_SI,          R_IR, N, X, Y, Beta, S0, I0, pi_IR, pi_S0, pi_I0,
                    Brate, delta, pinf_dogs, time, Itr, fclform)
    fc_1 <- I_SI_FC(I_SI_proposed, R_IR, N, X, Y, Beta, S0, I0, pi_IR, pi_S0, pi_I0,
                    Brate, delta, pinf_dogs, time, Itr, fclform)
    
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      I_SI <- I_SI_proposed
      AR_I_SI <- 1
    }
    else{
      I_SI <- I_SI
      AR_I_SI <- 0
    }
    return(list(I_SI, AR_I_SI))
  }
  
  
  Draw_R_IR <- function(R_IR, I_SI, N, Y, S0, I0, pi_IR, pi_I0, time, pro_prop){
    I <- calculate_I(I0, I_SI, R_IR)
    # if (any(I < 0)){return(list(-Inf, 0))}
    proposed_Trans <- Propose_New_Transition_Better_I_Promise_Sort_Of(
      R_IR, I, proposal_proportion = pro_prop)
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    R_IR_proposed <- proposed_Trans$compartment
    
    fc_0 <- R_IR_FC(R_IR,          I_SI, S0, I0, Y, N, pi_IR, pi_I0, time)
    fc_1 <- R_IR_FC(R_IR_proposed, I_SI, S0, I0, Y, N, pi_IR, pi_I0, time)
    checkNA(fc_0, "fc_0 in R_IR")
    checkNA(fc_1, "fc_1 in R_IR")
    
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      R_IR <- R_IR_proposed
      AR_R_IR <- 1
    }
    else{
      R_IR <- R_IR
      AR_R_IR <- 0
    }
    return(list(R_IR, AR_R_IR))
  }
  
  Draw_S0 <- function(S0, N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X, Brate, 
                      delta, pinf_dogs, Itr, time, fclform){
    proposed_Trans <- Propose_Initial_Compartment(N, S0, pi_S0)
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    S0_proposed <- proposed_Trans$Initial_Compartment
    
    fc_0 <- S0_FC(S0,          N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X, Brate, 
                  delta, pinf_dogs, Itr, time, fclform)
    fc_1 <- S0_FC(S0_proposed, N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X, Brate, 
                  delta, pinf_dogs, Itr, time, fclform)
    
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      S0 <- S0_proposed
      AR_S0 <- 1
    }
    else{
      S0 <- S0
      AR_S0 <- 0
    }
    return(list(S0, AR_S0))
  }
  
  Draw_I0 <- function(I0, S0, N, I_SI, R_IR, pi_I0, pi_IR, time){
    proposed_Trans <- Propose_Initial_Compartment(N, I0, pi_I0)
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    I0_proposed <- proposed_Trans$Initial_Compartment
    
    fc_0 <- I0_FC(I0,          S0, N, I_SI, R_IR, pi_I0, pi_IR, time)
    fc_1 <- I0_FC(I0_proposed, S0, N, I_SI, R_IR, pi_I0, pi_IR, time)
    
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      I0 <- I0_proposed
      AR_I0 <- 1
    }
    else{
      I0 <- I0
      AR_I0 <- 0
    }
    return(list(I0, AR_I0))
  }
  
  
  
  ## MCMC Algorithm:
  
  # Load data simulated from SIR model
  load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/SIR_generated_simulation.RData")
  # y.data <- sim_sir[sim_sir$Run==1,][,c(1,8:10)]
  # y.data <- y.data[y.data$Week!=0,]
  
  # Load data simulated from MLR model
  load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Simulation Study-SIR_vs_Multinomial/Multinomial_generated_simulation.RData")
  y.data <- sim_multinomial[sim_multinomial$Run==1,c(1:4)]
  week <- y.data$Week
  
  # Load epidemic process simulation 
  # load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/SIR_Simulation2.RData")
  data <- sim_sir[sim_sir$Run==1,][,c(1:7)]
  data <- data[data$Week!=0,]
  names(data) <- c("Week", names(data)[2:7])
  
  # Design matrix and starting betas for pi_SI function
  Itr <- 53
  X <- matrix(c(rep(1,length=Itr),
                c(dnorm(1:30, 15, 3)),rep(0, Itr-length(1:30))),ncol=2)
  # X <- matrix(c(rep(1, length=Itr),
  #             seq(1:Itr)), ncol=2)
  X <- X[1:max(y.data$Week),]
  
  # Truncate data at max week for which we have data to speed things up
  data <- data[1:nrow(X),]
  
  # Specify "known values"
  N <- 180000 # population size of Parnarim
  
  tp <- sum(0.510,0.019,0.350) #we're only including three of four categories - will need to re-weight by this
  pi_S0 <- 0.510/tp
  pi_I0 <- 0.019/tp
  pi_R0 <- 0.350/tp
  # pi_S0 <- 5/8
  # pi_I0 <- 1/8
  # pi_R0 <- 1/4
  
  # Initialize Betas
  # B <- matrix(c(0.1, -4))
  B <- matrix(c(0.25, -4))
  beta_mean <- c(0, 0)
  beta_sd <- c(0.2, 0.2)
  
  # Declare tuning parameters
  proposal_Mean <- 0
  proposal_SD <- 0.04
  pro_prop <- 0.20
  
  ## Indicate functional form of pi_SI
  fclform <- "exponential"
  delta <- 1 # For comparability, assuming only human contribution
  pinf_dogs <- 0.3 # This is not relevant for this analysis, but could be used in future work
  
  # We know that the rate of infection (Ab+, DTH- -- our tA group) is 26.1-27.4% for females
  # and 17.7-21.8% for males
  # Take extremes and say that 90% fall between 0.17 and 0.28 to determine prior parameter values
  source("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/derive_beta_parameters_from_quantiles.R")
  
  geom.mean6weeks <- qgeom(c(0.05, 0.5, 0.95), 0.10)
  geom.mean4weeks <- qgeom(c(0.05, 0.5, 0.95), 0.15)
  
  strong_prior_param_calc <- beta.parms.from.quantiles(q=c(0.10, 0.15), p=c(0.05, 0.95), plot=F)
  
  alpha_IR <- strong_prior_param_calc$a
  beta_IR <- strong_prior_param_calc$b
  
  # Set algorithm values
  n <- 1 # number of spatial locations
  TotalIterations <- 2000000
  Stride <- 100
  Tpts <- length(week)
  
  # Initial AR values
  AR_Beta <- 0
  AR_R_IR <- 0
  AR_I_SI <- 0
  AR_S0 <- 0
  AR_I0 <- 0
  
  # Set up storage for results
  Beta_Matrix <- matrix(data=NA, nrow=ncol(X), ncol=TotalIterations/Stride)
  
  pi_IR_Vector <- rep(NA, TotalIterations/Stride)
  
  I_SI_Array <- array(NA, dim=c(Tpts, n, TotalIterations/Stride))
  R_IR_Array <- array(NA, dim=c(Tpts, n, TotalIterations/Stride))
  
  pi_SI_Array <- array(NA, dim=c(Tpts, n, TotalIterations/Stride))
  
  S_Array <- array(NA, dim=c(Tpts, n, TotalIterations/Stride))
  I_Array <- array(NA, dim=c(Tpts, n, TotalIterations/Stride))
  R_Array <- array(NA, dim=c(Tpts, n, TotalIterations/Stride))
  
  S0_Vector <- rep(NA, TotalIterations/Stride)
  I0_Vector <- rep(NA, TotalIterations/Stride)
  
  # Initialize values
  Brate <- 30
  Beta <- matrix(B, ncol=ncol(X))
  pi_IR <- 0.10
  
  I_SI <- matrix(data$I_SI, ncol=n)
  R_IR <- matrix(data$R_IR, ncol=n)
  
  pi_SI <- matrix(data$pi_SI, ncol=n)
  
  S <- matrix(data$S, ncol=n)
  I <- matrix(data$I, ncol=n)
  R <- matrix(data$R, ncol=n)
  Y <- matrix(cbind(y.data$Y_S, y.data$Y_I, y.data$Y_R), ncol=3)
  
  # S0 <- floor(N*pi_S0)
  # I0 <- floor(N*pi_I0)
  S0 <- sim_sir[1,2]
  I0 <- sim_sir[1,3]
  R0 <- N-(S0+I0)
  
  ## Begin timing process
  ptm <- proc.time()
  
  for (iteration in 1:TotalIterations){
    
    # Draw Beta
    #   print("Beta")
    Beta_Info <- Draw_Beta(Beta=Beta, beta_mean=beta_mean, 
                           beta_sd=beta_sd, 
                           N=N, X=X, Y=Y, 
                           S0=S0, I0=I0,
                           I_SI=I_SI, R_IR=R_IR, pi_S0=pi_S0, pi_I0=pi_I0,
                           Brate=Brate, delta=delta, pinf_dogs=pinf_dogs, 
                           proposal_Mean=proposal_Mean, 
                           proposal_SD=proposal_SD, 
                           Itr=Itr, time=week, 
                           fclform=fclform)
    Beta <- Beta_Info[[1]]
    AR_Beta <- sum(c(AR_Beta, Beta_Info[[2]]))
    
    pi_SI <- calculate_pi_SI(X=X, Beta=Beta, I0=I0, I_SI=I_SI, R_IR=R_IR, 
                             N=N, Brate=Brate, delta=delta, 
                             pinf_dogs=pinf_dogs, Itr=Itr, fclform=fclform)
    
    # Draw I_SI
    #   print("I_SI")
    I_SI_Info <- Draw_I_SI(I_SI=I_SI, R_IR=R_IR, N=N, X=X, Y=Y, Beta=Beta, 
                           S0=S0, I0=I0, pi_IR=pi_IR, pi_S0=pi_S0, pi_I0=pi_I0,
                           Brate=Brate, delta=delta, pinf_dogs=pinf_dogs,
                           time=week, pro_prop=pro_prop, fclform=fclform)
    I_SI <- I_SI_Info[[1]]
    AR_I_SI <- sum(c(AR_I_SI, I_SI_Info[[2]]))
    S <- calculate_S(S0=S0, I_SI=I_SI)
    I <- calculate_I(I0=I0, I_SI=I_SI, R_IR=R_IR)
    R <- N-(S+I)
    
    pi_SI <- calculate_pi_SI(X=X, Beta=Beta, I0=I0, I_SI=I_SI, R_IR=R_IR,
                             N=N, Brate=Brate, delta=delta, 
                             pinf_dogs=pinf_dogs, Itr=Itr, fclform=fclform)
    
    # Draw R_IR
    #   print("R_IR")
    R_IR_Info <- Draw_R_IR(R_IR=R_IR, I_SI=I_SI, N=N, Y=Y, S0=S0, I0=I0,
                           pi_IR=pi_IR, pi_I0=pi_I0, time=week, pro_prop=pro_prop)
    R_IR <- R_IR_Info[[1]]
    AR_R_IR <- sum(c(AR_R_IR, R_IR_Info[[2]]))
    I <- calculate_I(I0=I0, I_SI=I_SI, R_IR=R_IR)
    R <- N-(S+I)
    
    #   print("Probs")
    pi_IR <- Draw_pi_IR(pi_IR=pi_IR, N=N, Y=Y, 
                        S0=S0, I0=I0, I_SI=I_SI, R_IR=R_IR, pi_I0=pi_I0,
                        alpha_IR=alpha_IR, beta_IR=beta_IR, time=week)
    
    # Draw New Initial Compartment Values
    S0_Info <- Draw_S0(S0=S0, N=N, I0=I0, I_SI=I_SI, R_IR=R_IR, pi_S0=pi_S0, pi_SI=pi_SI,
                       Beta=Beta, X=X, Brate=Brate, delta=delta, pinf_dogs=pinf_dogs, 
                       Itr=Itr, time=week, fclform=fclform)
    S0 <- S0_Info[[1]]
    AR_S0 <- sum(c(AR_S0, S0_Info[[2]]))
    
    I0_Info <- Draw_I0(I0=I0, S0=S0, N=N, I_SI=I_SI, R_IR=R_IR, pi_I0=pi_I0, pi_IR=pi_IR, time=week)
    I0 <- I0_Info[[1]]
    AR_I0 <- sum(c(AR_I0, I0_Info[[2]]))
    
    R0 <- N-(S0+I0)
    
    # Save the output somehow (start with beta), only want to save every Stride or 1000 iterations
    if (iteration %% Stride == 0)
    {
      cat(paste("Iteration: ", iteration, "\n", sep = ""))
      # output MCMC stuff here
      Beta_Matrix[,(iteration/Stride)] <- Beta
      
      pi_IR_Vector[(iteration/Stride)] <- pi_IR
      
      I_SI_Array[,,(iteration/Stride)] <- I_SI
      R_IR_Array[,,(iteration/Stride)] <- R_IR
      
      pi_SI_Array[,,(iteration/Stride)] <- pi_SI
      
      S_Array[,,(iteration/Stride)] <- S
      I_Array[,,(iteration/Stride)] <- I
      R_Array[,,(iteration/Stride)] <- R
      
      S0_Vector[(iteration/Stride)] <- S0
      I0_Vector[(iteration/Stride)] <- I0
    }
  }
  
  AR_Beta <- AR_Beta/TotalIterations
  AR_R_IR <- AR_R_IR/TotalIterations
  AR_I_SI <- AR_I_SI/TotalIterations
  AR_S0 <- AR_S0/TotalIterations
  AR_I0 <- AR_I0/TotalIterations
  
  proc.time() - ptm
  
  # fname <- paste0("SIR_simulation_fit", seed, "_Brate_", Brate, ".RData")
  # fname <- paste0("SIR_simulation_fit", seed, "_exp2", ".RData")
  # fname <- paste0("MLR_simulation_fit", seed, "_Brate_", Brate, ".RData")
  fname <- paste0("MLR_simulation_fit", seed, "_exp", ".RData")
  # fname <- paste0("Bayesian_final_proj_exp91417_", seed, ".RData")
  save(Beta_Matrix, AR_Beta, pi_IR_Vector, I_SI_Array, AR_I_SI, R_IR_Array, AR_R_IR,
       S0_Vector, AR_S0, I0_Vector, AR_I0, pi_SI_Array, S_Array, I_Array, R_Array, N, file=fname)
  return(TRUE)
}


## Run script:

runParallel <- TRUE

if (runParallel){
  library(parallel)
  Nchain <- 4
  
  cl <- makeCluster(Nchain)
  clusterExport(cl, "runMCMCChain")
  results <- parLapply(cl, 1:4, function(x){
    runMCMCChain(123123 + x)
  })
  stopCluster(cl)
} else{
  results <- lapply(1:4, function(x){
    runMCMCChain(123123 + x)
  })
}
