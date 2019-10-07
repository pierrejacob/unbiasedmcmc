#'@rdname sample_meetingtime
#'@title Sample coupled Markov chains until meeting
#'@description Sample two Markov chains, each following 'single_kernel' marginally,
#'until they meet, and report the meeting time, as well as the elapsed wall-clock time in seconds.
#'
#' This function does not record the trajectories of the chains, with the goal of being memory-light.
#' To record these trajectories, see \code{\link{sample_coupled_chains}}. To directly
#' compute unbiased estimators on the fly, see  \code{\link{sample_unbiasedestimator}}.
#'
#'@param single_kernel A list taking a state and returning a state, performing one step of a Markov kernel
#'@param coupled_kernel A list taking two states and returning two states, performing one step of a coupled Markov kernel;
#'it also returns a boolean "identical" indicating whether the two states are identical.
#'@param rinit A list representing the initial state of the chain, that can be given to 'single_kernel'
#'@param lag A time lag, equal to one by default
#'@param max_iterations A maximum number of iterations, at which to interrup the while loop; Inf by default
#'@return A list with
#'\itemize{
#'
#'\item meetingtime: the meeting time; equal to Inf if while loop was interrupted
#'
#'\item elapsedtime: elapsed wall-clock time, in seconds
#'}
#'@export
sample_meetingtime <- function(single_kernel, coupled_kernel, rinit, lag = 1, max_iterations = Inf){
  starttime <- Sys.time()
  # initialize
  state1 <- rinit(); state2 <- rinit()
  # move first chain
  time <- 0
  for (t in 1:lag){
    time <- time + 1
    state1 <- single_kernel(state1)
  }
  # move two chains until meeting (or until max_iterations)
  meetingtime <- Inf
  while (is.infinite(meetingtime) && (time < max_iterations)){
    time <- time + 1
    # use coupled kernel
    coupledstates <- coupled_kernel(state1, state2)
    state1 <- coupledstates$state1
    state2 <- coupledstates$state2
    # check if meeting happens
    if (coupledstates$identical) meetingtime <- time
  }
  currentime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(meetingtime = meetingtime, elapsedtime = elapsedtime))
}
