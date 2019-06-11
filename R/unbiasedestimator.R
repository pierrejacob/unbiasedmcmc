#'@rdname sample_unbiasedestimator
#'@title Unbiased MCMC estimators
#'@description Sample two Markov chains, each following 'single_kernel' marginally,
#' and 'coupled_kernel' jointly, until min(max(tau, m), max_iterations), where tau
#' is the first time the two chains meet (the "meeting time"). An unbiased estimator
#' of the expectation of a test function h is computed on the fly, between step k and step m, and returned.
#'
#' Allows for an arbitrary lag, i.e. X_t = Y_{t-lag}, and lag is one by default.
#'
#' Compared to \code{\link{sample_coupled_chains}} this function requires specifying
#' the test function, but does not record the trajectories, and thus is memory-light.
#'
#' If you're only interested in sampling meeting times, see \code{\link{sample_meetingtime}}.
#'
#'@param single_kernel A list taking a state and returning a state, performing one step of a Markov kernel
#'@param coupled_kernel A list taking two states and returning two states, performing one step of a coupled Markov kernel;
#'it also returns a boolean "identical" indicating whether the two states are identical.
#'@param rinit A list representing the initial state of the chain, that can be given to 'single_kernel'
#'@param h A test function of interest, which should take a chain state ("chain_state" entry of the output of "rinit", for instance)
#'and return a numeric vector
#'@param k An integer at which to start computing the unbiased estimator
#'@param m A time horizon: the chains are sampled until the maximum between m and the meeting time
#'@param lag A time lag, equal to one by default
#'@param max_iterations A maximum number of iterations, at which to interrup the while loop; Inf by default
#'@return A list with
#'\itemize{
#'
#'\item mcmcestimator: an MCMC estimator computed on the first chain, from step k to m
#'
#'\item correction: the bias correction term
#'
#'\item uestimator: unbiased estimator, equal to the sum of mcmcestimator and correction
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
sample_unbiasedestimator <- function(single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 1, lag = 1, max_iterations = Inf){
  if (k > m){
    print("error: k has to be less than m")
    return(NULL)
  }
  if (lag > m){
    print("error: lag has to be less than m")
    return(NULL)
  }
  starttime <- Sys.time()
  state1 <- rinit(); state2 <- rinit()
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(state1$chain_state)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(m-k+1, ceiling((t - k)/lag)) * (h(X_{t}) - h(Y_{t-lag})) for t=k+lag,..., tau - 1
  correction <- rep(0, dimh)
  time <- 0
  for (t in 1:lag){
    time <- time + 1
    state1 <- single_kernel(state1)
    if (time >= k){
      mcmcestimator <- mcmcestimator + h(state1$chain_state)
    }
  }
  if (time >= k + lag){
    correction <- correction + min(m-k+1, ceiling((time - k)/lag)) * (h(state1$chain_state) - h(state2$chain_state))
  }
  meetingtime <- Inf
  # time here is equal to lag; at this point we have X_lag,Y_0 and we are going to generate successively X_{t},Y_{t-lag} where time t is >= lag+1
  while ((time < max(meetingtime, m)) && (time < max_iterations)){
    time <- time + 1 # time is lag+1,lag+2,...
    if (is.finite(meetingtime)){
      state1 <- single_kernel(state1)
      state2 <- state1
      if (k <= time && time <= m){
        mcmcestimator <- mcmcestimator + h(state1$chain_state)
      }
    } else {
      res_coupled_kernel <- coupled_kernel(state1, state2)
      state1 <- res_coupled_kernel$state1
      state2 <- res_coupled_kernel$state2
      if (res_coupled_kernel$identical){
        meetingtime <- time
      }
      if (k <= time && time <= m){
        mcmcestimator <- mcmcestimator + h(state1$chain_state)
      }
      if (time >= k + lag){
        correction <- correction + min(m-k+1, ceiling((time - k)/lag)) * (h(state1$chain_state) - h(state2$chain_state))
      }
    }
  }
  uestimator <- mcmcestimator + correction
  cost <- lag + 2*(meetingtime - lag) + max(0, time - meetingtime)
  currentime <- Sys.time()
  elapsedtime <- as.numeric(as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(mcmcestimator = mcmcestimator / (m - k + 1), correction = correction / (m - k + 1), uestimator = uestimator / (m - k + 1),
              meetingtime = meetingtime, iteration = time, elapsedtime = elapsedtime, cost = cost))
}
