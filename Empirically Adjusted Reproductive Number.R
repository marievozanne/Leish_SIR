## Empirically adjusted reproduction number
## SIR compartmental model -- humans only

## This is the model fit with 2,000,000 total iterations and stride = 500
## File SIR multinomial_with_covariates_AND_sandfly_parallel_BPexp.R
## copy now included in Part 1 -- Real Data analysis folder
load("M:/CPHS/Projects/Leish/Marie/Marie Ozanne Code/Part 1 - Multinomial and SIR Comparison/Real_Data_Analysis-SIR_vs_Multinomial/MCMC_human_hist_w_cov_exp_123124.RData")

k <- dim(S_Array)[3]
n <- dim(S_Array)[1]

## Write function that calculates this:
R_EA_function <- function(S, I, pi_SI, pi_IR){
  if (is.null(dim(S))==FALSE) stop("Function valid for single location -- S must be a vector")
  if (is.null(dim(I))==FALSE) stop("Function valid for single location -- I must be a vector")
  exp_inf <- S*pi_SI
  avg_inf <- exp_inf/I
  G <- diag(avg_inf)
  n <- length(S)
  R_EA_comps <- G%*%((1-pi_IR)^c(0:(n-1)))
  R_EA_comps <- c(R_EA_comps,G[n,n]*(1-pi_IR)^(c(n:(n+20))))
  R_EA <- rep(NA, length(R_EA_comps))
  for (i in 1:length(R_EA_comps)){
    R_EA[i] <- sum(R_EA_comps[i:length(R_EA_comps)])
  }
  return(R_EA)
}

temp <- matrix(nrow=n+21, ncol=k)
for (i in 1:k){
  temp[,i] <- R_EA_function(S=S_Array[,,i], I=I_Array[,,i], 
                            pi_SI=pi_SI_Array[,,i], pi_IR=pi_IR_Vector[i])
}

## symmetric credible interval:
credible_intervals <- data.frame(t(apply(temp[,335:k], 1, quantile, c(0.025, 0.50, 0.975))))
credible_intervals$week <- c(0:48)

library(ggplot2)
library(reshape2)

m <- melt(credible_intervals, id.vars = "week", variable.name = "quantile")
m$bds <- rep(NA, nrow(m))
m[which(m$quantile!="X50."),]$bds <- "bounds"
m[which(m$quantile=="X50."),]$bds <- "median"

(p <- ggplot(credible_intervals[credible_intervals$week<29,], aes(x=week, y=X50.)) 
  + ylab(parse(text='R^(EA)')) + xlab("Week") 
  + scale_y_continuous(limits=c(0,3.4))
  + geom_ribbon(aes(ymin=X2.5., ymax=X97.5.), fill = "grey70")
  + geom_line(size=1.05) + theme_bw()
  + scale_x_continuous(breaks=seq(0, 28, 4))
  + geom_hline(yintercept = 1, lty="dashed", size=1.05, col="grey2")
  + theme(legend.position="none", axis.text=element_text(size=14), 
          axis.title = element_text(size=18)))

setEPS()
postscript("M:/CPHS/Projects/Leish/Marie/Papers in progress/Multinomial LR vs Compartmental Model -- Leish/Journal of Applied Statistics/Ozanne_Figure_4a.eps")
p
dev.off()

# tiff("M:/CPHS/Projects/Leish/Papers in progress/Multinomial LR vs Compartmental Model -- Leish/Ozanne_Figure_4.tiff", 
#      height = 6, width = 10, units = 'in', type="windows", res=800,
#      compression="lzw")
# p
# dev.off()