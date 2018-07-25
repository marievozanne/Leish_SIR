## Script for running Bayesian SIR model

runMCMCChain <- function(seed){
  set.seed(seed);
  setwd("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Real_Data_Analysis-SIR_vs_Multinomial/")
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
  Beta_FC <- function(Beta, beta_mean, beta_sd, N, X, Y, Gamma, Z, ind.time, 
                      S0, I0, I_SI, R_IR, pi_S0, pi_I0,
                      Brate, delta, pinf_dogs, Itr, time, fclform){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    
    proposed_pi_SI <- calculate_pi_SI(X, Beta, I0, I_SI, R_IR, N, Brate, 
                                      delta, pinf_dogs, Itr, fclform)
    
    if (any(proposed_pi_SI > 1)){ 
      return(-Inf)
    }
    
    # Individual level contributions
    ind.cov.contr <- data.frame(Z%*%Gamma, ind.time) 
    names(ind.cov.contr) <- c("Zgam", "Week")
    mean.Zgam <- ddply(ind.cov.contr, c("Week"), summarise, meanZgam <- mean(Zgam))
    names(mean.Zgam) <- c("Week", "meanZgam")
    
    I.aug <- I[time,]+mean.Zgam[,2]
    
    # Likelihood component:
    lik.data <- sum(dmultinomial(Y, prob = cbind(S[time,]/N, I.aug/N, 
                                                 1-(S[time,]+I.aug)/N), log=TRUE))
    lik.ic <- (dbinom(S0, N, pi_S0, log=TRUE) + 
                 dbinom(I0, N-S0, pi_I0/sum(c(pi_I0, 1-pi_S0-pi_I0)), log=TRUE))
    lik.tc <- sum(dbinom(I_SI, S, proposed_pi_SI, log = TRUE))
    prior <- sum(dnorm(Beta, mean = beta_mean, sd = beta_sd, log = TRUE))
    return(lik.data+lik.ic+lik.tc+prior)
  }
  
  ## Gamma_FC (to be used in incorporating individual-level covariates)
  Gamma_FC <- function(Gamma, gamma_mean, gamma_sd, Beta, N, X, Y, Z, S0, I0, I_SI, R_IR, 
                       pi_SI, pi_IR, pi_S0, pi_I0, Itr, time, ind.time){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    
    ## Individual level contributions
    ind.cov.contr <- data.frame(Z%*%Gamma, ind.time) 
    names(ind.cov.contr) <- c("Zgam", "Week")
    mean.Zgam <- ddply(ind.cov.contr, c("Week"), summarise, meanZgam <- mean(Zgam))
    names(mean.Zgam) <- c("Week", "meanZgam")
    
    I.aug <- I[time,]+mean.Zgam[,2]
    
    # Likelihood component:
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I.aug/N, 
                                               1-(S[time,]+I.aug)/N), log=TRUE))
    lik.ic <- (dbinom(S0, N, pi_S0, log=TRUE) + 
                 dbinom(I0, N, pi_I0, log=TRUE))
    lik.tc <- sum(dbinom(I_SI, S, pi_SI, log = TRUE)) +
      sum(dbinom(R_IR, I, pi_IR, log = TRUE))
    prior <- sum(dnorm(Gamma, mean = gamma_mean, sd = gamma_sd, log = TRUE))
    return(lik.data+lik.ic+lik.tc+prior)
  }
  
  ## pi_IR_FC
  pi_IR_FC <- function(pi_IR, S0, I0, I_SI, R_IR, N, Y, Z, Gamma, alpha_IR, beta_IR, time, ind.time){
    if (pi_IR < 0){
      return(-Inf)
    }
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    # if ((S0 + I0) > N){return(-Inf)}
    
    ## Individual level contributions
    ind.cov.contr <- data.frame(Z%*%Gamma, ind.time) 
    names(ind.cov.contr) <- c("Zgam", "Week")
    mean.Zgam <- ddply(ind.cov.contr, c("Week"), summarise, meanZgam <- mean(Zgam))
    names(mean.Zgam) <- c("Week", "meanZgam")
    
    I.aug <- I[time,]+mean.Zgam[,2]
    
    # Likelihood Component:
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I.aug/N, 
                                               1-(S[time,]+I.aug)/N), log = TRUE))
    lik.ic <- dbinom(I0, N, pi_I0, log=TRUE)
    lik.tc <- sum(dbinom(R_IR, I, pi_IR, log = TRUE))
    prior <- sum(dbeta(pi_IR, alpha_IR, beta_IR, log = TRUE))
    return(lik.data+lik.ic+lik.tc+prior)
  }
  
  ### --- Transition compartments --- 
  ## I_SI_FC
  I_SI_FC <- function(I_SI, R_IR, N, X, Y, Z, Beta, Gamma, S0, I0, pi_IR, pi_S0, pi_I0,
                      Brate, delta, pinf_dogs, time, ind.time, Itr, fclform){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    if (length(S0 + I0 > N) > 1){return("Uh-oh I_SI_FC")}
    # if ((S0 + I0) > N){return(-Inf)}
    
    pi_SI <- calculate_pi_SI(X, Beta, I0, I_SI, R_IR, N, Brate, delta, pinf_dogs, Itr, fclform)
    
    if (any(pi_SI > 1)){ ## Think about whether you really need this check
      return(-Inf)
    }
    
    ## Individual level contributions
    ind.cov.contr <- data.frame(Z%*%Gamma, ind.time) 
    names(ind.cov.contr) <- c("Zgam", "Week")
    mean.Zgam <- ddply(ind.cov.contr, c("Week"), summarise, meanZgam <- mean(Zgam))
    names(mean.Zgam) <- c("Week", "meanZgam")
    
    I.aug <- I[time,]+mean.Zgam[,2]
    
    #Likelihood component
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I.aug/N, 
                                               1-(S[time,]+I.aug)/N), log = TRUE))
    lik.ic <- (dbinom(S0, N, pi_S0, log=TRUE) + 
                 dbinom(I0, N-S0, pi_I0/sum(c(pi_I0,1-pi_S0-pi_I0)), log=TRUE))
    lik.tc <- sum(dbinom(I_SI, S, pi_SI, log = TRUE)) +
      sum(dbinom(R_IR, I, pi_IR, log = TRUE))
    return(lik.data+lik.ic+lik.tc)
  }
  
  ## R_IR_FC 
  R_IR_FC <- function(R_IR, I_SI, S0, I0, Y, Z, Gamma, N, pi_IR, pi_I0, time, ind.time){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    # if ((S0 + I0) > N){return(-Inf)}
    
    ## Individual level contributions
    ind.cov.contr <- data.frame(Z%*%Gamma, ind.time) 
    names(ind.cov.contr) <- c("Zgam", "Week")
    mean.Zgam <- ddply(ind.cov.contr, c("Week"), summarise, meanZgam <- mean(Zgam))
    names(mean.Zgam) <- c("Week", "meanZgam")
    
    I.aug <- I[time,]+mean.Zgam[,2]
    
    #Likelihood component
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I.aug/N, 
                                               1-(S[time,]+I.aug)/N), log = TRUE))
    lik.ic <- dbinom(I0, N-S0, pi_I0/sum(c(pi_I0,1-pi_S0-pi_I0)), log=TRUE)
    lik.tc <- sum(dbinom(R_IR, I, pi_IR, log = TRUE))
    return(lik.data+lik.ic+lik.tc)
  }
  
  ### --- Initial compartments ---
  S0_FC <- function(S0, N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X, Gamma, Z,
                    Brate, delta, pinf_dogs, Itr, time, fclform){
    S <- calculate_S(S0, I_SI)
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(S < 0)){return(-Inf)}
    if (any(I < 0)){return(-Inf)}
    
    pi_SI <- calculate_pi_SI(X, Beta, I0, I_SI, R_IR, N, 
                             Brate, delta, pinf_dogs, Itr, fclform)
    if (any(pi_SI > 1)){return(-Inf)}
    
    # Individual level contributions
    ind.cov.contr <- data.frame(Z%*%Gamma, ind.time) 
    names(ind.cov.contr) <- c("Zgam", "Week")
    mean.Zgam <- ddply(ind.cov.contr, c("Week"), summarise, meanZgam <- mean(Zgam))
    names(mean.Zgam) <- c("Week", "meanZgam")
    
    I.aug <- I[time,]+mean.Zgam[,2]
    
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I.aug/N, 
                                               1-(S[time,]+I.aug)/N), log = TRUE))
    lik.ic <- sum(dbinom(S0, N, pi_S0, log=TRUE))
    lik.tc <- sum(dbinom(I_SI, S, pi_SI, log=TRUE))
    return(lik.data+lik.ic+lik.tc)
  }
  
  I0_FC <- function(I0, S0, N, I_SI, R_IR, pi_I0, pi_S0, pi_IR, Gamma, Z, time){
    I <- calculate_I(I0, I_SI, R_IR)
    
    if (any(I < 0)){return(-Inf)}
    
    # Individual level contributions
    ind.cov.contr <- data.frame(Z%*%Gamma, ind.time) 
    names(ind.cov.contr) <- c("Zgam", "Week")
    mean.Zgam <- ddply(ind.cov.contr, c("Week"), summarise, meanZgam <- mean(Zgam))
    names(mean.Zgam) <- c("Week", "meanZgam")
    
    I.aug <- I[time,]+mean.Zgam[,2]
    
    lik.data <- sum(dmultinomial(Y, prob=cbind(S[time,]/N, I.aug/N, 
                                               1-(S[time,]+I.aug)/N), log = TRUE))
    lik.ic <- sum(dbinom(I0, N-S0, pi_I0/sum(c(pi_I0,1-pi_S0-pi_I0)), log=TRUE))
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
    
    proposalProb <- ifelse(From_Compartment[idx,] == 0, 
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
  Draw_Beta <- function(Beta, beta_mean, beta_sd, Gamma, N, X, Y, Z, 
                        S0, I0, I_SI, R_IR, pi_S0, pi_I0,
                        Brate, delta, pinf_dogs, proposal_Mean, proposal_SD, 
                        Itr, time, ind.time, fclform){
    Beta_proposed <- Propose_Coef(Beta, proposal_Mean, proposal_SD)
    
    # Beta, beta_mean, beta_sd, N, X, Y, Gamma, Z, ind.time, 
    # S0, I0, I_SI, R_IR, pi_S0, pi_I0,
    # Brate, delta, pinf_dogs, Itr, time, fclform
    
    fc_0 <- Beta_FC(Beta,          beta_mean, beta_sd, N, X, Y, Gamma, Z, ind.time,
                    S0, I0, I_SI, R_IR, pi_S0, pi_I0,
                    Brate, delta, pinf_dogs, Itr, time, fclform)
    fc_1 <- Beta_FC(Beta_proposed, beta_mean, beta_sd, N, X, Y, Gamma, Z, ind.time,
                    S0, I0, I_SI, R_IR, pi_S0, pi_I0, 
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
  
  Draw_Gamma <- function(Gamma, gamma_mean, gamma_sd, Beta, N, X, Y, Z, S0, I0, 
                         I_SI, R_IR, pi_SI, pi_IR, pi_I0, pi_S0,
                         proposal_Mean, proposal_SD, Itr, time, ind.time){
    ## Propose new gamma values
    Gamma_proposed <- Propose_Coef(Gamma, proposal_Mean, proposal_SD)
    ## Full conditionals
    fc_0 <- Gamma_FC(Gamma,          gamma_mean, gamma_sd, Beta, N, X, Y, Z, 
                     S0, I0, I_SI, R_IR, 
                     pi_SI, pi_IR, pi_S0, pi_I0, Itr, time, ind.time)
    fc_1 <- Gamma_FC(Gamma_proposed, gamma_mean, gamma_sd, Beta, N, X, Y, Z, 
                     S0, I0, I_SI, R_IR, 
                     pi_SI, pi_IR, pi_S0, pi_I0, Itr, time, ind.time)
    ## MH Step
    Rat <- fc_1 - fc_0
    if (log(runif(1)) < min(Rat, 0) && !is.nan(Rat)){
      Gamma <- Gamma_proposed
      AR_Gamma <- 1
    }
    else{
      Gamma <- Gamma
      AR_Gamma <- 0
    }
    ## Return parameter value and indicator of whether accepted
    return(list(Gamma, AR_Gamma))
  }
  
  Draw_pi_IR <- function(pi_IR, N, Y, S0, I0, I_SI, R_IR, pi_I0, alpha_IR, beta_IR, time){
    I <- calculate_I(I0, I_SI, R_IR)
    if (any(I < 0)){return(-Inf)}
    ## Draw new value of pi_IR from posterior -- Gibbs step
    pi_IR_proposed <- Propose_New_Transition_Probability(R_IR, I, alpha_IR, beta_IR)
    pi_IR_proposed
  }
  
  Draw_I_SI <- function(I_SI, R_IR, N, X, Y, Z, Beta, Gamma, S0, I0, pi_IR, pi_S0, pi_I0,
                        Brate, delta, pinf_dogs, time, ind.time, pro_prop, fclform){
    ## Calculate compartment size
    S <- calculate_S(S0, I_SI)
    
    ## Propose new transition compartment
    proposed_Trans <- Propose_New_Transition_Better_I_Promise_Sort_Of(
      I_SI, S, proposal_proportion = pro_prop)
    if (any(proposed_Trans$Compartment > S)){
      print("Propose_New_Transition_Better_I_Promise_Sort_Of is not doing its job.")
    }
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    I_SI_proposed <- proposed_Trans$compartment
    
    ## Full conditionals
    fc_0 <- I_SI_FC(I_SI,          R_IR, N, X, Y, Z, Beta, Gamma, S0, I0, pi_IR, pi_S0, pi_I0,
                    Brate, delta, pinf_dogs, time, ind.time, Itr, fclform)
    fc_1 <- I_SI_FC(I_SI_proposed, R_IR, N, X, Y, Z, Beta, Gamma, S0, I0, pi_IR, pi_S0, pi_I0,
                    Brate, delta, pinf_dogs, time, ind.time, Itr, fclform)
    
    ## MH Step
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      I_SI <- I_SI_proposed
      AR_I_SI <- 1
    }
    else{
      I_SI <- I_SI
      AR_I_SI <- 0
    }
    ## Return transition compartment and indicator of whether accepted
    return(list(I_SI, AR_I_SI))
  }
  
  Draw_R_IR <- function(R_IR, I_SI, N, Y, Z, Gamma, S0, I0, pi_IR, pi_I0, 
                        time, ind.time, pro_prop){
    ## Calculate compartment size
    I <- calculate_I(I0, I_SI, R_IR)
    
    ## Propose new transition compartment
    proposed_Trans <- Propose_New_Transition_Better_I_Promise_Sort_Of(
      R_IR, I, proposal_proportion = pro_prop)
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    R_IR_proposed <- proposed_Trans$compartment
    
    ## Full conditionals
    fc_0 <- R_IR_FC(R_IR,          I_SI, S0, I0, Y, Z, Gamma, N, pi_IR, pi_I0, time, ind.time)
    fc_1 <- R_IR_FC(R_IR_proposed, I_SI, S0, I0, Y, Z, Gamma, N, pi_IR, pi_I0, time, ind.time)
    checkNA(fc_0, "fc_0 in R_IR")
    checkNA(fc_1, "fc_1 in R_IR")
    
    ## MH Step
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      R_IR <- R_IR_proposed
      AR_R_IR <- 1
    }
    else{
      R_IR <- R_IR
      AR_R_IR <- 0
    }
    ## Return transition compartment and indicator of whether accepted
    return(list(R_IR, AR_R_IR))
  }
  
  Draw_S0 <- function(S0, N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X,  Gamma, Z,
                      Brate, delta, pinf_dogs, Itr, time, fclform){
    ## Propose new initial compartment
    proposed_Trans <- Propose_Initial_Compartment(N, S0, pi_S0)
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    S0_proposed <- proposed_Trans$Initial_Compartment
    
    ## Full conditionals
    fc_0 <- S0_FC(S0,          N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X, Gamma, Z,
                  Brate, delta, pinf_dogs, Itr, time, fclform)
    fc_1 <- S0_FC(S0_proposed, N, I0, I_SI, R_IR, pi_S0, pi_SI, Beta, X, Gamma, Z, 
                  Brate, delta, pinf_dogs, Itr, time, fclform)
    
    ## MH Step
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      S0 <- S0_proposed
      AR_S0 <- 1
    }
    else{
      S0 <- S0
      AR_S0 <- 0
    }
    ## Return initial compartment and indicator of whether accepted
    return(list(S0, AR_S0))
  }
  
  Draw_I0 <- function(I0, S0, N, I_SI, R_IR, pi_I0, pi_S0, pi_IR, Gamma, Z, time){
    ## Propose new initial compartment 
    proposed_Trans <- Propose_Initial_Compartment(N, I0, pi_I0)
    
    p_prp_gvn_cur <- proposed_Trans$prob_proposed
    p_cur_gvn_prp <- proposed_Trans$prob_current
    I0_proposed <- proposed_Trans$Initial_Compartment
    
    ## Full conditionals
    fc_0 <- I0_FC(I0,          S0, N, I_SI, R_IR, pi_I0, pi_S0, pi_IR, Gamma, Z, time)
    fc_1 <- I0_FC(I0_proposed, S0, N, I_SI, R_IR, pi_I0, pi_S0, pi_IR, Gamma, Z, time)
    
    ## MH Step
    Rat <- ((fc_1 + p_cur_gvn_prp) - (fc_0 + p_prp_gvn_cur))
    if (log(runif(1)) < min(Rat,0) && is.na(Rat)==FALSE){
      I0 <- I0_proposed
      AR_I0 <- 1
    }
    else{
      I0 <- I0
      AR_I0 <- 0
    }
    ## Return initial compartment and indicator of whether accepted
    return(list(I0, AR_I0))
  }
  
  ## MCMC Algorithm:
  
  # Load observed data
  load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/human_cons_no_cov.RData")
  load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/human_new_2.RData")
  y.data <- human_cons_no_cov
  week <- y.data$Week
  
  # Load epidemic process simulation
  load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/SIR_Simulation2.RData")
  data <- sim_mat2[sim_mat2$Run==2,]
  names(data) <- c("Week", names(data)[2:12])
  
  # Design matrix and starting betas for pi_SI function
  Itr <- 53
  X <- matrix(c(rep(1,length=Itr),
                c(dnorm(1:30, 15, 3)),rep(0, Itr-length(1:30))),ncol=2)
  # X <- matrix(c(rep(1, length=Itr),
  #               seq(1:Itr)), ncol=2)
  X <- X[1:max(y.data$Week),]
  
  # Initialize betas
  B <- matrix(c(1,-4))
  beta_mean <- c(0,-1)
  beta_sd <- c(0.001,0.001)
  
  # Truncate data at max week for which we have data to speed things up
  data <- data[1:nrow(X),]
  
  # Design matrix and starting gammas for individual covariates component
  Z <- model.matrix(Disease_category ~ Sex + Cat_Age + Area + Week, human_new)
  ind.time <- Z[,ncol(Z)]
  Z <- Z[,-c(1,ncol(Z))]
  
  
  # Initialize gammas
  Gamma <- matrix(c(1,1.2,0.8,-0.6)+rnorm(4, 0, 0.1), nrow=ncol(Z)) 
  gamma_mean <- c(0,0,0,0)
  gamma_sd <- c(0.008,0.008,0.008,0.008)
  
  # Declare tuning parameters
  N <- 180000
  
  tp <- sum(0.510,0.019,0.350) #we're only including three of four categories - will need to re-weight by this
  pi_S0 <- 0.510/tp
  pi_I0 <- 0.019/tp
  
  S0 <- floor(N*pi_S0)
  I0 <- floor(N*pi_I0)
  R0 <- N-(S0+I0)
  
  proposal_Mean_Beta <- 0
  proposal_SD_Beta <- 0.003
  proposal_Mean_Gamma <- 0
  proposal_SD_Gamma <- 0.01
  
  pro_prop <- 0.28
  
  ## Indicate functional form of pi_SI
  fclform <- "exponential"
  delta <- 1 # For comparability, assuming only human contribution
  pinf_dogs <- 0.3 # This is not relevant for this analysis, but could be used in future work
  
  # We know that the rate of infection (Ab+, DTH- -- our tA group) is 26.1-27.4% for females
  # and 17.7-21.8% for males
  # Take extremes and say that 90% fall between 0.17 and 0.28 to determine prior parameter values
  source("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/derive_beta_parameters_from_quantiles.R")
  
  # p1=0.10 ## median of 6 weeks
  # p2=0.15 ## median of 4 weeks
  
  # qgeom(c(0.05, 0.5, 0.95), 0.10)
  # qgeom(c(0.05, 0.5, 0.95), 0.15)
  
  # geom.mean13weeks <- qgeom(c(0.05, 0.5, 0.95), 0.05)
  # geom.mean03weeks <- qgeom(c(0.05, 0.5, 0.95), 0.20)
  
  # p1=0.05 ## median of 13 weeks
  # p2=0.20 ## median of 03 weeks
  
  qgeom(c(0.05, 0.5, 0.95), 0.07) ## median of 9 weeks
  qgeom(c(0.05, 0.5, 0.95), 0.20) ## median of 3 weeks
  
  p1=0.07 ## median of 9 weeks
  p2=0.20 ## median of 3 weeks
  
  strong_prior_param_calc <- beta.parms.from.quantiles(q=c(p1, p2), p=c(0.15, 0.85), plot=F)
  
  alpha_IR <- strong_prior_param_calc$a
  beta_IR <- strong_prior_param_calc$b
  
  # alpha_IR <- 1
  # beta_IR <- 1
  
  # Set algorithm values
  n <- 1 # number of spatial locations
  TotalIterations <- 2000000
  # TotalIterations <- 2000
  Stride <- 500
  Tpts <- nrow(data)
  
  # Initial AR values
  AR_Beta <- 0
  AR_Gamma <- 0
  AR_R_IR <- 0
  AR_I_SI <- 0
  AR_S0 <- 0
  AR_I0 <- 0
  
  # Set up storage for results
  Beta_Matrix <- matrix(data=NA, nrow=ncol(X), ncol=TotalIterations/Stride)
  Gamma_Matrix <- matrix(data=NA, nrow=ncol(Z), ncol=TotalIterations/Stride)
  
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
  Brate <- 20
  Beta <- matrix(B, ncol=ncol(X))
  pi_IR <- 0.1098373 
  
  I_SI <- matrix(data$I_SI, ncol=n)
  R_IR <- matrix(data$R_IR, ncol=n)
  
  pi_SI <- matrix(data$pi_SI, ncol=n)
  
  S <- matrix(data$S, ncol=n)
  I <- matrix(data$I, ncol=n)
  R <- matrix(data$R, ncol=n)
  Y <- matrix(cbind(y.data$S, y.data$I, y.data$R), ncol=3)
  
  ptm <- proc.time()
  
  for (iteration in 1:TotalIterations){
    # Draw Initial Compartment Values
    S0_Info <- Draw_S0(S0=S0, N=N, I0=I0, I_SI=I_SI, R_IR=R_IR, pi_S0=pi_S0, pi_SI=pi_SI,
                       Beta=Beta, X=X, Gamma=Gamma, Z=Z,
                       Brate=Brate, delta=delta, pinf_dogs=pinf_dogs, 
                       Itr=Itr, time=week, fclform=fclform)
    S0 <- S0_Info[[1]]
    AR_S0 <- sum(c(AR_S0, S0_Info[[2]]))
    
    I0_Info <- Draw_I0(I0=I0, S0=S0, N=N, I_SI=I_SI, R_IR=R_IR, 
                       pi_I0=pi_I0, pi_S0=pi_S0, pi_IR=pi_IR, 
                       Gamma=Gamma, Z=Z, time=week)
    I0 <- I0_Info[[1]]
    AR_I0 <- sum(c(AR_I0, I0_Info[[2]]))
    
    R0 <- N-(S0+I0)
    
    # Draw Beta
    #   print("Beta")
    Beta_Info <- Draw_Beta(Beta=Beta, beta_mean=beta_mean, 
                           beta_sd=beta_sd, Gamma=Gamma, 
                           N=N, X=X, Y=Y, 
                           Z=Z, S0=S0, I0=I0,
                           I_SI=I_SI, R_IR=R_IR, pi_S0=pi_S0, pi_I0=pi_I0,
                           Brate=Brate, delta=delta, pinf_dogs=pinf_dogs, 
                           proposal_Mean=proposal_Mean_Beta, 
                           proposal_SD=proposal_SD_Beta, 
                           Itr=Itr, time=week, 
                           ind.time=ind.time, fclform=fclform)
    Beta <- Beta_Info[[1]]
    AR_Beta <- sum(c(AR_Beta, Beta_Info[[2]]))
    
    pi_SI <- calculate_pi_SI(X=X, Beta=Beta, I0=I0, I_SI=I_SI, R_IR=R_IR, 
                             N=N, Brate=Brate, delta=delta, 
                             pinf_dogs=pinf_dogs, Itr=Itr, fclform=fclform)
    
    # Draw Gamma 
    Gamma_Info <- Draw_Gamma(Gamma=Gamma, gamma_mean=gamma_mean, 
                             gamma_sd=gamma_sd, Beta=Beta, 
                             N=N, X=X, Y=Y, Z=Z, 
                             S0=S0, I0=I0, I_SI=I_SI, R_IR=R_IR,
                             pi_SI=pi_SI, pi_IR=pi_IR, pi_S0=pi_S0, pi_I0=pi_I0,
                             proposal_Mean = proposal_Mean_Gamma, 
                             proposal_SD = proposal_SD_Gamma, 
                             Itr = Itr, time=week, ind.time=ind.time)
    Gamma <- Gamma_Info[[1]]
    AR_Gamma <- sum(c(AR_Gamma, Gamma_Info[[2]]))
    
    # Draw I_SI
    #   print("I_SI")
    I_SI_Info <- Draw_I_SI(I_SI=I_SI, R_IR=R_IR, N=N, X=X, Y=Y, Z=Z, Beta=Beta, Gamma=Gamma,
                           S0=S0, I0=I0, pi_IR=pi_IR, pi_S0=pi_S0, pi_I0=pi_I0,
                           Brate=Brate, delta=delta, pinf_dogs=pinf_dogs,
                           time=week, ind.time=ind.time, pro_prop=pro_prop, fclform=fclform)
    I_SI <- I_SI_Info[[1]]
    AR_I_SI <- sum(c(AR_I_SI, I_SI_Info[[2]]))
    
    S <- calculate_S(S0=S0, I_SI=I_SI)
    I <- calculate_I(I0=I0, I_SI=I_SI, R_IR=R_IR)
    R <- N-(S+I)
    pi_SI <- calculate_pi_SI(X=X, Beta=Beta, I0=I0, I_SI=I_SI, R_IR=R_IR, N=N, Brate=Brate, 
                             delta=delta, pinf_dogs=pinf_dogs, Itr=Itr, fclform=fclform)
    
    # Draw R_IR
    #   print("R_IR")
    R_IR_Info <- Draw_R_IR(R_IR=R_IR, I_SI=I_SI, N=N, Y=Y, Z=Z, Gamma=Gamma, S0=S0, I0=I0, 
                           pi_IR=pi_IR, pi_I0=pi_I0, time=week, ind.time=ind.time, pro_prop=pro_prop)
    R_IR <- R_IR_Info[[1]]
    AR_R_IR <- sum(c(AR_R_IR, R_IR_Info[[2]]))
    
    I <- calculate_I(I0=I0, I_SI=I_SI, R_IR=R_IR)
    R <- N-(S+I)
    
    #   print("Probs")
    pi_IR <- Draw_pi_IR(pi_IR=pi_IR, N=N, Y=Y, 
                        S0=S0, I0=I0, I_SI=I_SI, R_IR=R_IR, 
                        alpha_IR=alpha_IR, beta_IR=beta_IR, time=week)
    
    # Save the output somehow (start with beta), only want to save every Stride or 1000 iterations
    if (iteration %% Stride == 0)
    {
      cat(paste("Iteration: ", iteration, "\n", sep = ""))
      # output MCMC stuff here
      Beta_Matrix[,(iteration/Stride)] <- Beta
      Gamma_Matrix[,(iteration/Stride)] <- Gamma
      
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
  AR_Gamma <- AR_Gamma/TotalIterations
  AR_R_IR <- AR_R_IR/TotalIterations
  AR_I_SI <- AR_I_SI/TotalIterations
  AR_S0 <- AR_S0/TotalIterations
  AR_I0 <- AR_I0/TotalIterations
  
  proc.time() - ptm

  # fname <- paste0("MCMC_human_hist_w_cov3_exp_", seed, ".RData")
  fname <- paste0("MCMC_human_hist_w_cov_weakp_1585_", seed, ".RData")
  save(Beta_Matrix, AR_Beta, Gamma_Matrix, AR_Gamma, pi_IR_Vector, 
       I_SI_Array, AR_I_SI,
       R_IR_Array, AR_R_IR, pi_SI_Array, S_Array, I_Array, R_Array, 
       I0_Vector, S0_Vector, N, AR_S0, AR_I0, file=fname)
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
