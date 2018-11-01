### This script plays with coupled parallel tempering and unbiased estimators
### for a basic Ising model, with different values of the temperatures
### i.e. coupled Gibbs sampler at each temperature and coupled swap moves between chains
library(debiasedmcmc)
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
registerDoParallel(cores = 2)

ising_unbiased_estimator <- function(betas, proba_swapmove, k = 0, m = 1, max_iterations = Inf){
  nchains  <- length(betas)
  ss_ <- c(-4,-2,0,2,4)
  probas_ <-  sapply(betas, function(beta) exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta)))
  # initialize
  chain_states1 <- ising_pt_rinit(nchains)
  chain_states2 <- ising_pt_rinit(nchains)
  sumstates1 <- unlist(lapply(chain_states1, debiasedmcmc:::ising_sum_))
  sumstates2 <- unlist(lapply(chain_states2, debiasedmcmc:::ising_sum_))
  # mcmcestimator computes the natural statistic for each chain
  mcmcestimator <- sumstates1
  if (k > 0){
    mcmcestimator <- rep(0, nchains)
  }
  # move first chain
  iter <- 1
  res_single_kernel <- ising_pt_single_kernel(chain_states1, sumstates1, betas, probas_, proba_swapmove)
  chain_states1 <- res_single_kernel$chain_states
  sumstates1 <- res_single_kernel$sumstates
  # correction term computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, nchains)
  if (k == 0){
    correction <- correction + min(1, (0 - k + 1)/(m - k + 1) )  * (sumstates1 - sumstates2)
  }
  # accumulate mcmc estimator
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + sumstates1
  }
  # iterate
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      # only need to use single kernel after meeting
      res_ <- ising_pt_single_kernel(chain_states1, sumstates1, betas, probas_, proba_swapmove)
      chain_states1 <- res_$chain_states
      sumstates1 <- res_$sumstates
      chain_states2 <- chain_states1
      sumstates2 <- sumstates1
      # accumulate mcmc estimator
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + sumstates1
      }

    } else {
      # use coupled kernel
      res_ <- ising_pt_coupled_kernel(chain_states1, chain_states2, sumstates1, sumstates2, betas, probas_, proba_swapmove)
      chain_states1 <- res_$chain_states1
      chain_states2 <- res_$chain_states2
      sumstates1 <- res_$sumstates1
      sumstates2 <- res_$sumstates2
      # check if meeting happens
      allequal <- all(sapply(1:nchains, function(i) all(chain_states1[[i]] == chain_states2[[i]])))
      if (allequal && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
      if (k <= iter){
        # accumulate mcmc estimator
        if (iter <= m){
          mcmcestimator <- mcmcestimator + sumstates1
        }
        # accumulate correction term
        correction <- correction + min(1, (iter-1 - k + 1)/(m - k + 1) ) * (sumstates1 - sumstates2)
      }
    }
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  # compute unbiased estimator
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction
  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}


proba_swapmove <- 0.01
nchains <- 10
betas <- seq(from = 0.30, to = 0.40, length.out = nchains)

# nrep <- 50
# res_preliminary <- foreach(i = 1:nrep) %dorng% {
#   ising_unbiased_estimator(betas, proba_swapmove, k = 0, m = 1)
# }
# hist(sapply(res_preliminary, function(x) x$meetingtime))
# floor(as.numeric(quantile(sapply(res_preliminary, function(x) x$meetingtime), probs = 0.95)))

k <- 300
nrep <- 250
#
# results <- foreach(i = 1:nrep) %dorng% {
#   ising_unbiased_estimator(betas, proba_swapmove, k = k, m = 5*k)
# }
# save(k, nrep, betas, nchains, proba_swapmove, results, file = "ising.testunbiasedestimators.RData")
load("ising.testunbiasedestimators.RData")
nrep <- length(results)

meetings <- sapply(results, function(x) x$meetingtime)
hist(meetings)
k
m <- 5*k

uestimators <- t(sapply(results, function(x) x$uestimator))
m <- colMeans(uestimators)
s <- apply(uestimators, 2, sd) / sqrt(nrep)
df <- data.frame(betas = betas, m = m, s = s)
g <- ggplot(df, aes(x = betas, y = m)) + geom_line() + geom_point()
g <- g + geom_errorbar(aes(ymin = m - 2*s, ymax = m + 2*s))
g

### Compare with usual MCMC
niterations <- 15000
# history
history_sumstates <- matrix(0, ncol = nchains, nrow = niterations)
# precomputed probability for single-site flips
ss_ <- c(-4,-2,0,2,4)
probas_ <-  sapply(betas, function(beta) exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta)))
# initialization
current_states <- ising_pt_rinit(nchains)
sumstates <- unlist(lapply(current_states, debiasedmcmc:::ising_sum_))
history_sumstates[1,] <- sumstates
nswap_attempts <- 0
nswap_accepts <- rep(0, nchains-1)
# MCMC loop
for (iteration in 2:niterations){
  res_ <- ising_pt_single_kernel(current_states, sumstates, betas, probas_, proba_swapmove)
  current_states <- res_$chain_states
  sumstates <- res_$sumstates
  nswap_accepts <- nswap_accepts + res_$nswap_accepts
  nswap_attempts <- nswap_attempts + res_$nswap_attempts
  history_sumstates[iteration,] <- sumstates
}

# swap acceptance in %
cat(paste0(round(100*nswap_accepts/nswap_attempts, 2), " %"), "\n")
par(mfrow = c(2,1))
# traceplots
matplot(apply(history_sumstates[2000:niterations,1:5], 2, function(x) cumsum(x)/(1:(niterations-2000+1))), type = "l")
# plot estimated average natural parameters
plot(betas, colMeans(history_sumstates[2000:niterations,]), type = "l")

par(mfrow = c(1,1))
matplot(history_sumstates[2000:niterations,5:10], type = "l")
library(coda)
effectiveSize(history_sumstates[5000:niterations,])
spectrum0.ar(history_sumstates[5000:niterations,])
##

mcmc_estimates <- colMeans(history_sumstates[5000:niterations,])
df$mcmc <- mcmc_estimates

g <- ggplot(df, aes(x = betas, y = m)) + geom_line() + geom_point()
g <- g + geom_errorbar(aes(ymin = m - 2*s, ymax = m + 2*s))
g + geom_line(aes(y = mcmc), colour = "red", linetype = 2)

