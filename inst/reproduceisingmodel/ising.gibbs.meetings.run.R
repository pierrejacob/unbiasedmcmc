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
#
# res_ <- debiasedmcmc:::coupled_gibbs_sweep(chain_state1, chain_state2, proba_)
# all(res_$state1 == res_$state2)

for (iter in 2:niterations){
  res_ <- ising_coupled_kernel(chain_state1, chain_state2, proba_)
  chain_state1 <- res_$chain_state1
  chain_state2 <- res_$chain_state2
  sumstates1[iter] <- debiasedmcmc:::ising_sum_(chain_state1)
  sumstates2[iter] <- debiasedmcmc:::ising_sum_(chain_state2)
}
matplot(cbind(sumstates1, sumstates2), type = "l")

coupled_chains <- function(beta, m = 1, max_iterations = Inf){
  ss_ <- c(-4,-2,0,2,4)
  proba_ <-  exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta))
  chain_state1 <- ising_rinit()
  chain_state2 <- ising_rinit()
  # dimstate <- dim(chain_state1)[1]
  # sumstates1 <- rep(0, m+preallocate+1)
  # sumstates2 <- rep(0, m+preallocate)
  # samples1 <- array(dim = c(m+preallocate+1, dimstate, dimstate))
  # samples2 <- array(dim = c(m+preallocate, dimstate, dimstate))
  # nrowsamples1 <- m+preallocate+1
  # sumstates1[1] <- debiasedmcmc:::ising_sum_(chain_state1)
  # sumstates2[1] <- debiasedmcmc:::ising_sum_(chain_state2)
  current_nsamples1 <- 1
  chain_state1 <- ising_single_kernel(chain_state1, proba_)
  current_nsamples1 <- current_nsamples1 + 1
  # sumstates1[current_nsamples1] <- debiasedmcmc:::ising_sum_(chain_state1)
  # samples1[current_nsamples1,,] <- chain_state1
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
    # if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      # new_rows <- nrowsamples1-1
      # nrowsamples1 <- nrowsamples1 + new_rows
      # sumstates1 <- c(sumstates1, rep(0, new_rows))
      # sumstates2 <- c(sumstates2, rep(0, new_rows))
    # }
    # sumstates1[current_nsamples1+1] <- debiasedmcmc:::ising_sum_(chain_state1)
    # sumstates2[current_nsamples1] <- debiasedmcmc:::ising_sum_(chain_state2)
    # samples1[current_nsamples1+1,,] <- chain_state1
    # samples2[current_nsamples1,,] <- chain_state2
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  # sumstates1 <- sumstates1[1:current_nsamples1]
  # sumstates2 <- sumstates2[1:(current_nsamples1-1)]
  # samples1 <- samples1[1:current_nsamples1,,,drop=F]
  # samples2 <- samples2[1:(current_nsamples1-1),,,drop=F]
  return(list(meetingtime = meetingtime, iteration = iter, finished = finished))
}


# nrep <- 100
# beta1 <- 0.2
# ccs1_prelim <- foreach(irep = 1:nrep) %dorng% {
#   coupled_chains(beta1, m = 1, max_iterations = Inf, preallocate = 10)
# }
# sapply(ccs1_prelim, function(x) x$meetingtime) %>% summary
#

nrep <- 100
betas <- seq(from = 0.3, to = 0.55, length.out = 15)
singlesite.meetings.df <- data.frame()
for (ibeta in seq_along(betas)){
  print(ibeta)
  beta <- betas[ibeta]
  ccs_ <- foreach(irep = 1:nrep) %dorng% {
    coupled_chains(beta, m = 1, max_iterations = Inf)
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
