#'@rdname dinversegamma
#'@title compute log-density of inverse gamma
#'@description at x > 0 and with given parameters alpha, beta, given by
#'  alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'@export
dinversegamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

#'@rdname rinversegamma
#'@title Sample from inverse gamma
#'@description with given parameters alpha, beta, with log-density given by
#' alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'@export
rinversegamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}

#'@rdname rinversegamma_coupled
#'@title Sample from maximally coupled inverse gamma
#'@description with given parameters alpha1, alpha2, beta1, beta2,
#'where the parametrization is that the log-density of IG(alpha, beta) is
#' alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'@export
rinversegamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  f <- get_max_coupling(function(n) rinversegamma(n, alpha1, beta1),
                        function(x) dinversegamma(x, alpha1, beta1),
                        function(n) rinversegamma(n, alpha2, beta2),
                        function(x) dinversegamma(x, alpha2, beta2))
  return(f())
}
