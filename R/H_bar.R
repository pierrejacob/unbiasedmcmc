#'@rdname H_bar
#'@title Compute unbiased estimators from coupled chains
#'@description Compute the proposed unbiased estimators, for each of the element
#'in the list 'c_chains'. The integral of interest is that of the function h,
#'which can be multivariate. The estimator uses the variance reduction technique
#'whereby the estimator is the MCMC average between times k and m, with probability
#'going to one as k increases.
#'@export
H_bar <- function(c_chains, h = function(x) x, k = 0, m = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if (m > maxiter){
    print("error: m has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # infer the lag from the number of rows in samples1 and samples2
  lag <- dim(c_chains$samples1)[1] - dim(c_chains$samples2)[1]
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  h_of_chain <- apply(X = c_chains$samples1[(k+1):(m+1),,drop=F], MARGIN = 1, FUN = h)
  if (is.null(dim(h_of_chain))){
    h_of_chain <- matrix(h_of_chain, ncol = 1)
  } else {
    h_of_chain <- t(h_of_chain)
  }
  H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
  # next, add bias correction terms
  # Delta_t refers to h(X_t) - h(Y_{t-lag})
  if (c_chains$meetingtime <= k + lag){
    # nothing else to add, because Delta_t = 0 for all t >= meeting time
  } else {
    for (time in (k+lag):(c_chains$meetingtime-1)){
      # time is the index t of X_{t} where the chain start from X_{0}
      coefficient_t <- min(m-k+1, ceiling((time-k)/lag))
      Delta_t <- h(c_chains$samples1[time+1,]) - h(c_chains$samples2[time-lag+1,])
      H_bar <- H_bar + coefficient_t * Delta_t
    }
  }
  # return result divided by number of terms i.e. m - k + 1
  return(H_bar / (m - k + 1))
}
