# Run coupled chains until max(tau, m) where tau is the meeting time and m specified by user
#'@rdname sample_coupled_chains
#'@title Sample coupled Markov chains
#'@description sample two Markov chains, each following 'single_kernel' marginally,
#' and 'coupled_kernel' jointly, until min(max(tau, m), max_iterations), where tau
#' is the first time the two chains meet. Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}.
#'@export
sample_coupled_chains <- function(single_kernel, coupled_kernel, rinit, m = 1, lag = 1, max_iterations = Inf, preallocate = 10){
  starttime <- Sys.time()
  state1 <- rinit(); state2 <- rinit()
  dimstate <- length(state1$chain_state)
  nrowsamples1 <- m+preallocate+lag
  samples1 <- matrix(nrow = nrowsamples1, ncol = dimstate)
  samples2 <- matrix(nrow = nrowsamples1-lag, ncol = dimstate)
  samples1[1,] <- state1$chain_state
  samples2[1,] <- state2$chain_state
  # current_nsamples1 <- 1
  time <- 0
  for (t in 1:lag){
    time <- time + 1
    state1 <- single_kernel(state1)
    samples1[time+1,] <- state1$chain_state
  }
  # current_nsamples1 <- current_nsamples1 + 1
  # iter <- 1
  meetingtime <- Inf
  while ((time < max(meetingtime, m)) && (time < max_iterations)){
    time <- time + 1 # time is lag+1,lag+2,...
    if (is.finite(meetingtime)){
      state1 <- single_kernel(state1)
      state2 <- state1
    } else {
      res_coupled_kernel <- coupled_kernel(state1, state2)
      state1 <- res_coupled_kernel$state1
      state2 <- res_coupled_kernel$state2
      if (res_coupled_kernel$identical){
        meetingtime <- time
      }
    }
    if ((time+1) > nrowsamples1){
      new_rows <- nrowsamples1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = dimstate))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = dimstate))
    }
    samples1[time+1,] <- state1$chain_state
    samples2[time-lag+1,] <-   state2$chain_state
  }
  samples1 <- samples1[1:(time+1),,drop=F]
  samples2 <- samples2[1:(time-lag+1),,drop=F]
  cost <- lag + 2*(meetingtime - lag) + max(0, time - meetingtime)
  currentime <- Sys.time()
  elapsedtime <- as.numeric(as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(samples1 = samples1, samples2 = samples2,
              meetingtime = meetingtime, iteration = time, elapsedtime = elapsedtime, cost = cost))
}

#' ## function to continue coupled chains until step m
#' ## c_chain should be the output of coupled_chains
#' ## and m should be more than c_chain$iteration, otherwise returns c_chain
#' #'@rdname continue_coupled_chains
#' #'@title Continue coupled MCMC chains up to m steps
#' #'@description ## function to continue coupled chains until step m
#' #' c_chain should be the output of coupled_chains
#' #' and m should be more than c_chain$iteration, otherwise returns c_chain
#' #'@export
#' continue_coupled_chains <- function(c_chain, single_kernel, m = 1, ...){
#'   if (m <= c_chain$iteration){
#'     ## nothing to do
#'     return(c_chain)
#'   } else {
#'     niterations <- m - c_chain$iteration
#'     chain_state1 <- c_chain$samples1[c_chain$iteration+1,]
#'     dimstate <- length(chain_state1)
#'     samples1 <- matrix(nrow = niterations, ncol = dimstate)
#'     samples2 <- matrix(nrow = niterations, ncol = dimstate)
#'     for (iteration in 1:niterations){
#'       chain_state1 <- single_kernel(chain_state1, ...)
#'       samples1[iteration,] <- chain_state1
#'       samples2[iteration,] <- chain_state1
#'     }
#'     c_chain$samples1 <- rbind(c_chain$samples1, samples1)
#'     c_chain$samples2 <- rbind(c_chain$samples2, samples2)
#'     c_chain$iteration <- m
#'     return(c_chain)
#'   }
#' }

