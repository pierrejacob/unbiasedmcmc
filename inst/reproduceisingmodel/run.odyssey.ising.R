## Example of file to be run in batch mode on the cluster
## This one runs parallel tempering for the ising model,
## produces unbiased estimators until "TIME" is up
## with a prescribed k and m
## The test function is the natural statistic of the lattice,
## and no history of chains is recorded (so memory requirement is low)

# time allocated to this
TIME <- 10*3600
# TIME <- 12

library(debiasedmcmc)

ising_unbiased_estimator <- function(betas, proba_swapmove, k = 0, m = 1, max_iterations = Inf, totalduration = Inf){
  ptm <- proc.time()
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
  # Check if time is up already
  elapsedtime <- as.numeric((proc.time() - ptm)[3])
  if (elapsedtime > totalduration){
    # time's up
    return(list(finished = FALSE, message = "interrupted because time's up"))
  }
  # iterate
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    elapsedtime <- as.numeric((proc.time() - ptm)[3])
    if (elapsedtime > totalduration){
      # time's up
      return(list(finished = FALSE, message = "interrupted because time's up"))
    }
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


library(parallel)
# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

# inverse temperatures of interest
nchains <- 16
betas <- seq(from = 0.3, to = 0.55, length.out = nchains)
proba_swapmove <- 0.01
k <- 50000
m <- 250000

ptm_irep <- proc.time()
last_start <- ptm_irep
nsamples_ <- 0
durations_ <- c()
elapsedtime <- 0
timeleft <- TIME
results_ <- list()
while(elapsedtime < TIME){
  if (nsamples_ == 0){ # then produce a sample
    ## TO CHANGE
    res <- try(ising_unbiased_estimator(betas, proba_swapmove, k = k, m = m, max_iterations = Inf))
    if (inherits(res, "try-error")) res <- list(finished = FALSE)
  } else { # then produce a sample if time permits
    ## TO CHANGE
    res <- try(ising_unbiased_estimator(betas, proba_swapmove, k = k, m = m, max_iterations = Inf, totalduration = timeleft))
    if (inherits(res, "try-error")) res <- list(finished = FALSE)
  }
  elapsedtime <- as.numeric((proc.time() - ptm_irep)[3])
  timeleft <- TIME - elapsedtime
  duration_ <- as.numeric((proc.time() - last_start)[3])
  last_start <- proc.time()
  if (res$finished){
    nsamples_ <- nsamples_ + 1
    results_[[nsamples_]] <- res
    durations_ <- c(durations_, duration_)
  }
  save(results_, nsamples_, durations_, file = paste0("output/ising.N", nchains, ".k", k, ".m", m, ".ID", igrid, ".RData"))
}

# durations_
# nsamples_
# sapply(results_, function(x) x$meetingtime)
