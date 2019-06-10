## sample meeting time
## arguments:
# single_kernel takes an object and returns an object
# coupled_kernel takes two objects and returns two objects "state1" and "state2", as well as an indicator "identical" (=0 or 1) indicating whether the two states are identical
# rinit takes no argument and returns an object, to be used as an input for "single_kernel" or "coupled_kernel"
# lag specifies the desired lag between the two chains
# max_iterations specifies when to stop the while loop
## returns:
# the meeting time (time at which the two chains coincide) and the elapsed wall-clock time in seconds
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
  elapsedtime <- as.numeric(as.duration(lubridate::ymd_hms(currentime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(meetingtime = meetingtime, elapsedtime = elapsedtime))
}
