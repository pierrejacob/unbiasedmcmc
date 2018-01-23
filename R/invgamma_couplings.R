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
  f <- get_max_coupling(function(n) rigamma(n, alpha1, beta1),
                        function(x) digamma(x, alpha1, beta1),
                        function(n) rigamma(n, alpha2, beta2),
                        function(x) digamma(x, alpha2, beta2))
  return(f())
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

