
# This file consists of the general-purpose functions coupled_chains, continue_coupled_chains,
# and H_bar, which implement our debiased mcmc algorithm for general kernels and functions h(.)

# from coupled_chains -----------------------------------------------------
# Run coupled chains until max(tau, K) where tau is the meeting time and K specified by user
#'@export
coupled_chains <- function(single_kernel, coupled_kernel, rinit, ..., K = 1, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  p <- length(chain_state1)
  samples1 <- matrix(nrow = K+11, ncol = p)
  samples2 <- matrix(nrow = K+10, ncol = p)
  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  current_nsamples1 <- 1
  chain_state1 <- single_kernel(chain_state1, ...)
  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      chain_state1 <- single_kernel(chain_state1, ...)
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, ...)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrow(samples1)){
      print('increase nrow')
      new_rows <- nrow(samples2)
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = ncol(samples1)))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = ncol(samples2)))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}


## function to continue coupled chains until step K
## c_chain should be the output of coupled_chains
## and K should be more than c_chain$iteration, otherwise returns c_chain
#'@export
continue_coupled_chains <- function(c_chain, single_kernel, K = 1, ...){
  if (K <= c_chain$iteration){
    ## nothing to do
    return(c_chain)
  } else {
    niterations <- K - c_chain$iteration
    chain_state1 <- c_chain$samples1[c_chain$iteration+1,]
    p <- length(chain_state1)
    samples1 <- matrix(nrow = niterations, ncol = p)
    samples2 <- matrix(nrow = niterations, ncol = p)
    for (iteration in 1:niterations){
      chain_state1 <- single_kernel(chain_state1, ...)
      samples1[iteration,] <- chain_state1
      samples2[iteration,] <- chain_state1
    }
    c_chain$samples1 <- rbind(c_chain$samples1, samples1)
    c_chain$samples2 <- rbind(c_chain$samples2, samples2)
    c_chain$iteration <- K
    return(c_chain)
  }
}



# from h_bar --------------------------------------------------------------

#'@export
H_bar <- function(c_chains, h = function(x) x, k = 0, K = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if (K > maxiter){
    print("error: K has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  h_of_chain <- apply(X = c_chains$samples1[(k+1):(K+1),,drop=F], MARGIN = 1, FUN = h)
  if (is.null(dim(h_of_chain))){
    h_of_chain <- matrix(h_of_chain, ncol = 1)
  } else {
    h_of_chain <- t(h_of_chain)
  }
  H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
  if (c_chains$meetingtime <= k + 1){
    # nothing else to add
  } else {
    deltas <- matrix(0, nrow = maxiter - k + 1, ncol = p)
    deltas_term <- rep(0, p)
    for (t in k:min(maxiter-1, c_chains$meetingtime-1)){ # t is as in the report, where the chains start at t=0
      coefficient <- min(t - k + 1, K - k + 1)
      delta_tp1 <- h(c_chains$samples1[t + 1 + 1,]) - h(c_chains$samples2[t+1,]) # the +1's are because R starts indexing at 1
      deltas_term <- deltas_term + coefficient * delta_tp1
    }
    H_bar <- H_bar + deltas_term
  }
  return(H_bar / (K - k + 1))
}

# Note that we can find 'H_bar_old' in the original h_bar.R file in inst

