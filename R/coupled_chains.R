# Run coupled chains until max(tau, K) where tau is the meeting time and K specified by user
#'@rdname coupled_chains
#'@title Coupled MCMC chains
#'@description sample two MCMC chains, each following 'single_kernel' marginally,
#' and 'coupled_kernel' jointly, until min(max(tau, m), max_iterations), where tau
#' is the first time the two chains meet. Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}.
#'@export
coupled_chains <- function(single_kernel, coupled_kernel, rinit, ..., m = 1, max_iterations = Inf, preallocate = 10){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  dimstate <- length(chain_state1)
  samples1 <- matrix(nrow = m+preallocate+1, ncol = dimstate)
  samples2 <- matrix(nrow = m+preallocate, ncol = dimstate)
  nrowsamples1 <- m+preallocate+1
  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  current_nsamples1 <- 1
  sres1 <- single_kernel(chain_state1, ...)
  chain_state1 <- sres1$state
  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      sres1 <- single_kernel(chain_state1, ...)
      chain_state1 <- sres1$state
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, ...)
      chain_state1 <- res_coupled_kernel$state1
      chain_state2 <- res_coupled_kernel$state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1-1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = dimstate))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = dimstate))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

#'@export
unbiasedestimator <- function(single_kernel, coupled_kernel, rinit, ..., h = function(x) x, k = 0, m = 1, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(chain_state1)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  sres1 <- single_kernel(chain_state1, ...)
  chain_state1 <- sres1$state
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(chain_state1)
  }
  # current_nsamples1 <- current_nsamples1 + 1
  # samples1[current_nsamples1,] <- chain_state1
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      sres1 <- single_kernel(chain_state1, ...)
      chain_state1 <- sres1$state
      chain_state2 <- chain_state1
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(chain_state1)
      }
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, ...)
      chain_state1 <- res_coupled_kernel$state1
      chain_state2 <- res_coupled_kernel$state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
      }
    }
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction
  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

## function to continue coupled chains until step m
## c_chain should be the output of coupled_chains
## and m should be more than c_chain$iteration, otherwise returns c_chain
#'@rdname continue_coupled_chains
#'@title Continue coupled MCMC chains up to m steps
#'@description ## function to continue coupled chains until step m
#' c_chain should be the output of coupled_chains
#' and m should be more than c_chain$iteration, otherwise returns c_chain
#'@export
continue_coupled_chains <- function(c_chain, single_kernel, m = 1, ...){
  if (m <= c_chain$iteration){
    ## nothing to do
    return(c_chain)
  } else {
    niterations <- m - c_chain$iteration
    chain_state1 <- c_chain$samples1[c_chain$iteration+1,]
    dimstate <- length(chain_state1)
    samples1 <- matrix(nrow = niterations, ncol = dimstate)
    samples2 <- matrix(nrow = niterations, ncol = dimstate)
    for (iteration in 1:niterations){
      chain_state1 <- single_kernel(chain_state1, ...)
      samples1[iteration,] <- chain_state1
      samples2[iteration,] <- chain_state1
    }
    c_chain$samples1 <- rbind(c_chain$samples1, samples1)
    c_chain$samples2 <- rbind(c_chain$samples2, samples2)
    c_chain$iteration <- m
    return(c_chain)
  }
}

