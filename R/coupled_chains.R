#'@rdname sample_coupled_chains
#'@title Sample coupled Markov chains
#'@description Sample two Markov chains, each following 'single_kernel' marginally,
#' and 'coupled_kernel' jointly, until min(max(tau, m), max_iterations), where tau
#' is the first time the two chains meet (the "meeting time").
#'
#' Or more precisely, they meet with a delay of lag, i.e. X_t = Y_{t-lag}, and lag is one by default.
#'
#' Once the coupled chains are obtained, unbiased estimators can be computed for arbitrary test
#' functions via the function \code{\link{H_bar}}.
#'
#' If you're only interested in sampling meeting times, see \code{\link{sample_meetingtime}}.
#'
#'
#'@param single_kernel A list taking a state and returning a state, performing one step of a Markov kernel
#'@param coupled_kernel A list taking two states and returning two states, performing one step of a coupled Markov kernel;
#'it also returns a boolean "identical" indicating whether the two states are identical.
#'@param rinit A list representing the initial state of the chain, that can be given to 'single_kernel'
#'@param m A time horizon: the chains are sampled until the maximum between m and the meeting time
#'@param lag A time lag, equal to one by default
#'@param max_iterations A maximum number of iterations, at which to interrup the while loop; Inf by default
#'@param preallocate A number of anticipated iterations, used to pre-allocate memory; 10 by default
#'@return A list with
#'\itemize{
#'
#'\item samples1: the first chain, of length max(m, tau)
#'
#'\item samples2: the second chain, of length max(m, tau) - lag
#'
#'\item meetingtime: the meeting time; equal to Inf if while loop was interrupted
#'
#'\item iteration: final iteration; could be equal to m, to meetingtime, or to max_iterations
#'
#'\item elapsedtime: elapsed wall-clock time, in seconds
#'
#'\item cost: computing cost in terms of calls to Markov kernels (counting coupled kernel as twice the cost)
#'}
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
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(samples1 = samples1, samples2 = samples2,
              meetingtime = meetingtime, iteration = time, elapsedtime = elapsedtime, cost = cost))
}

