
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
