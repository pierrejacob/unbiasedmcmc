#'@rdname digamma
#'@title compute log-density of inverse gamma
#'@description at x > 0 and with given parameters alpha, beta, given by
#'  alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'@export
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

#'@rdname rigamma
#'@title Sample from inverse gamma
#'@description with given parameters alpha, beta, with log-density given by
#' alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'@export
rigamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}

#'@rdname rigamma_coupled
#'@title Sample from maximally coupled inverse gamma
#'@description with given parameters alpha1, alpha2, beta1, beta2,
#'where the parametrization is that the log-density of IG(alpha, beta) is
#' alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'@export
rigamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  x <- rigamma(1, alpha1, beta1)
  if (digamma(x, alpha1, beta1) + log(runif(1)) < digamma(x, alpha2, beta2)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rigamma(1, alpha2, beta2)
      reject <- (digamma(y, alpha2, beta2) + log(runif(1)) < digamma(y, alpha1, beta1))
    }
    return(c(x,y))
  }
}

#' #'@rdname rigamma_transport_coupled
#' #'@title Sample from transport coupled inverse gamma
#' #'@description with given parameters alpha, beta1, beta2
#' #'@export
#' rigamma_transport_coupled <- function(alpha, beta1, beta2){
#'   u <- runif(1)
#'   Z1 <- 1/qgamma(p = u, shape = alpha, rate = beta1)
#'   Z2 <- 1/qgamma(p = u, shape = alpha, rate = beta2)
#'   return(list(Z1 = Z1, Z2 = Z2))
#' }

