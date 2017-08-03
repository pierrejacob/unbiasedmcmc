#'@rdname digamma
#'@title compute log-density of inverse gamma
#'@description at x > 0 and with given parameters alpha, beta
#'@export
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

#'@rdname rigamma
#'@title Sample from inverse gamma
#'@description with given parameters alpha, beta
#'@export
rigamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}

#'@rdname rigamma_coupled
#'@title Sample from coupled inverse gamma
#'@description with given parameters alpha, beta1, beta2
#'@export
rigamma_coupled <- function(alpha, beta1, beta2){
  # note: check below is unnecessary
  if (beta1 < beta2){
    Z <- rigamma(1, alpha, beta1)
    logratio <- -(beta2 - beta1) / Z
    if (log(runif(1)) < logratio){
      return(list(Z1 = Z, Z2 = Z))
    } else {
      return(list(Z1 = Z, Z2 = rigamma(1, alpha, beta2)))
    }
  } else {
    Z <- rigamma(1, alpha, beta2)
    logratio <- -(beta1 - beta2) / Z
    if (log(runif(1)) < logratio){
      return(list(Z1 = Z, Z2 = Z))
    } else {
      return(list(Z1 = rigamma(1, alpha, beta1), Z2 = Z))
    }
  }
}

#'@rdname rigamma_transport_coupled
#'@title Sample from transport coupled inverse gamma
#'@description with given parameters alpha, beta1, beta2
#'@export
rigamma_transport_coupled <- function(alpha, beta1, beta2){
  u <- runif(1)
  Z1 <- 1/qgamma(p = u, shape = alpha, rate = beta1)
  Z2 <- 1/qgamma(p = u, shape = alpha, rate = beta2)
  return(list(Z1 = Z1, Z2 = Z2))
}

