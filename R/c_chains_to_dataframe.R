#'@export
c_chains_to_dataframe <- function(c_chains, k, m){
  nsamples <- length(c_chains)
  approximation <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    measure <- unbiasedmcmc:::get_measure_(c_chains[[irep]], k, m)
    data.frame(rep = rep(irep, length(measure$weights)), weights = measure$weights, atoms = measure$atoms)
  }
  approximation$weights <- approximation$weights / nsamples
  approximation <- data.frame(unbiasedmcmc:::prune_(as.matrix(approximation[do.call(order, as.list(approximation[,3:ncol(approximation)])),])))
  approximation <- approximation[approximation$weights != 0,]
  return(approximation)
}
