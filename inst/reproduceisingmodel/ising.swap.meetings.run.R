### This script plays with coupled parallel tempering
### for a basic Ising model, with different values of the temperatures
### i.e. coupled Gibbs sampler at each temperature and coupled swap moves between chains
library(debiasedmcmc)
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
registerDoParallel(cores = detectCores()-2)

# size of the grid
size <- 32

ising_pt_coupled_chains <- function(betas, proba_swapmove, m = 1, max_iterations = Inf, preallocate = 10){
  nchains  <- length(betas)
  nswap_attempts <- 0
  nswap_accepts1 <- rep(0, nchains-1)
  ss_ <- c(-4,-2,0,2,4)
  probas_ <-  sapply(betas, function(beta) exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta)))
  # initialization
  chain_states1 <- ising_pt_rinit(nchains)
  chain_states2 <- ising_pt_rinit(nchains)
  sumstates1 <- unlist(lapply(chain_states1, debiasedmcmc:::ising_sum_))
  sumstates2 <- unlist(lapply(chain_states2, debiasedmcmc:::ising_sum_))
  history_sumstates1 <- matrix(0, ncol = nchains, nrow = m+preallocate+1)
  history_sumstates2 <- matrix(0, ncol = nchains, nrow = m+preallocate)
  nrowsamples1 <- m+preallocate+1
  history_sumstates1[1,] <- sumstates1
  history_sumstates2[1,] <- sumstates2
  current_nsamples1 <- 1
  # one step of first chain
  res_ <- ising_pt_single_kernel(chain_states1, sumstates1, betas, probas_, proba_swapmove)
  nswap_attempts <- nswap_attempts + res_$nswap_attempts
  nswap_accepts1 <- nswap_accepts1 + res_$nswap_accepts
  chain_states1 <- res_$chain_states
  sumstates1 <- res_$sumstates
  current_nsamples1 <- current_nsamples1 + 1
  history_sumstates1[current_nsamples1,] <- sumstates1
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){ # advance first chain and copy
      res_ <- ising_pt_single_kernel(chain_states1, sumstates1, betas, probas_, proba_swapmove)
      nswap_attempts <- nswap_attempts + res_$nswap_attempts
      nswap_accepts1 <- nswap_accepts1 + res_$nswap_accepts
      chain_states1 <- res_$chain_states
      sumstates1 <- res_$sumstates
      chain_states2 <- chain_states1
      sumstates2 <- sumstates1
    } else { # advance coupled chains
      res_ <- ising_pt_coupled_kernel(chain_states1, chain_states2, sumstates1, sumstates2, betas, probas_, proba_swapmove)
      nswap_attempts <- nswap_attempts + res_$nswap_attempts
      nswap_accepts1 <- nswap_accepts1 + res_$nswap_accepts1
      chain_states1 <- res_$chain_states1
      chain_states2 <- res_$chain_states2
      sumstates1 <- res_$sumstates1
      sumstates2 <- res_$sumstates2
      allequal <- all(sapply(1:nchains, function(i) all(chain_states1[[i]] == chain_states2[[i]])))
      if (allequal && !meet){ # record meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      new_rows <- nrowsamples1-1
      nrowsamples1 <- nrowsamples1 + new_rows
      history_sumstates1 <- rbind(history_sumstates1, matrix(0, nrow = new_rows, ncol = nchains))
      history_sumstates2 <- rbind(history_sumstates2, matrix(0, nrow = new_rows, ncol = nchains))
    }
    history_sumstates1[current_nsamples1+1,] <- sumstates1
    history_sumstates2[current_nsamples1,]   <- sumstates2
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  history_sumstates1 <- history_sumstates1[1:current_nsamples1,]
  history_sumstates2 <- history_sumstates2[1:(current_nsamples1-1),]
  return(list(history_sumstates1 = history_sumstates1, history_sumstates2 = history_sumstates2,
              meetingtime = meetingtime, iteration = iter, finished = finished,
              nswap_attempts = nswap_attempts, nswap_accepts1 = nswap_accepts1))
}


## proba_swapmove <- 0.01
## nchains <- 10
## betas <- seq(from = 0.35, to = 0.45, length.out = nchains)
## res_ <- ising_pt_coupled_chains(betas, proba_swapmove)
## res_$meetingtime
## cat(100*res_$nswap_accepts1 / res_$nswap_attempts, "%\n")

nrep <- 50
nchains_values <- c(8, 12, 16)
# nchains_values <- c(4,8,12,16,24,32)
proba_swapmove <- 0.01
for (nchains in nchains_values){
  print(nchains)
  betas <- seq(from = 0.3, to = 0.55, length.out = nchains)
  ccs_ <- foreach(irep = 1:nrep) %dorng% {
    ising_pt_coupled_chains(betas, proba_swapmove)
  }
  meetings <- sapply(ccs_, function(x) x$meetingtime)
  save(nchains, nrep, betas, meetings, proba_swapmove, file = paste0("ising.swapN", nchains, ".meetings.RData"))
}


nchains.df <- data.frame()
for (nchains in nchains_values){
  load(paste0("ising.swapN", nchains, ".meetings.RData"))
  nchains.df <- rbind(nchains.df, data.frame(nchains = nchains, mean = mean(meetings),
                                             q90 = as.numeric(quantile(meetings, probs = 0.9)),
                                             max = max(meetings)))
}
nchains.df

#
setmytheme()
ggplot(nchains.df, aes(x = nchains, y = mean)) + geom_line() + geom_point() + scale_y_log10(breaks = c(1e3,1e4,1e5), limits = c(1e4, 1e6)) +
  scale_x_continuous(breaks = nchains_values)

