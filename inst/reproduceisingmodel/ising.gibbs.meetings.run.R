### This script plays with a coupled Gibbs sampler (i.e. single site updates)
### for a basic Ising model, with different values of the temperatures
library(debiasedmcmc)
rm(list = ls())
set.seed(21)
#
library(dplyr)
library(doRNG)
library(doParallel)
registerDoParallel(cores = detectCores()-2)

library(abind)

# size of the grid
size <- 32
# possible values of sum of neighbors
ss_ <- c(-4,-2,0,2,4)
# inverse temperature
beta <- 0.4
# precomputed probability for single-site flips
proba_ <-  exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta))
# initialization
chain_state1 <- ising_rinit()
chain_state2 <- ising_rinit()
# MCMC loop
niterations <- 1000
# store history of sum of states
sumstates1 <- rep(0, niterations)
sumstates1[1] <- debiasedmcmc:::ising_sum_(chain_state1)
sumstates2 <- rep(0, niterations)
sumstates2[1] <- debiasedmcmc:::ising_sum_(chain_state2)

for (iter in 2:niterations){
  res_ <- ising_coupled_kernel(chain_state1, chain_state2, proba_)
  chain_state1 <- res_$chain_state1
  chain_state2 <- res_$chain_state2
  sumstates1[iter] <- debiasedmcmc:::ising_sum_(chain_state1)
  sumstates2[iter] <- debiasedmcmc:::ising_sum_(chain_state2)
}

matplot(cbind(sumstates1, sumstates2), type = "l")

ising_meeting <- function(beta, m = 1, max_iterations = Inf){
  ss_ <- c(-4,-2,0,2,4)
  proba_ <-  exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta))
  chain_state1 <- ising_rinit()
  chain_state2 <- ising_rinit()
  current_nsamples1 <- 1
  chain_state1 <- ising_single_kernel(chain_state1, proba_)
  current_nsamples1 <- current_nsamples1 + 1
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      chain_state1 <- ising_single_kernel(chain_state1, proba_)
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- ising_coupled_kernel(chain_state1, chain_state2, proba_)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  return(list(meetingtime = meetingtime, iteration = iter, finished = finished))
}


nrep <- 10
# betas <- seq(from = 0.3, to = 0.55, length.out = 15)
betas <- c(0.3, 0.317857142857143, 0.335714285714286, 0.353571428571429,
  0.371428571428571, 0.389285714285714, 0.407142857142857, 0.425,
  0.442857142857143, 0.460714285714286, 0.478571428571429)
singlesite.meetings.df <- data.frame()
for (ibeta in seq_along(betas)){
  print(ibeta)
  beta <- betas[ibeta]
  ccs_ <- foreach(irep = 1:nrep) %dorng% {
    ising_meeting(beta, m = 1, max_iterations = Inf)
  }
  meetingtimes <- sapply(ccs_, function(x) x$meetingtime)
  singlesite.meetings.df <- rbind(singlesite.meetings.df,
                                  data.frame(meeting = meetingtimes, beta = rep(beta, nrep)))
  print(summary(meetingtimes))
  save(nrep, betas, singlesite.meetings.df, file = "ising.singlesite.meetings.RData")
}

load(file = "ising.singlesite.meetings.RData")
tail(singlesite.meetings.df)

setmytheme()
g <- ggplot(singlesite.meetings.df %>% group_by(beta) %>% summarise(m = mean(meeting)),
       aes(x = beta, y = m)) + geom_line() + geom_point() + scale_y_log10()
g <- g + xlab(expression(theta)) + ylab("average meeting time")
g

# ggsave(filename = "ising.singlesite.meetings.pdf", plot = g, width = 8, height = 6)
