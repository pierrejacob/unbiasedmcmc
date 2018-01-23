#'@rdname rgamma_coupled
#'@title Sample from maximally coupled Gamma
#'@description with given parameters alpha1, alpha2, beta1, beta2,
#'where the parametrization is that beta is the rate, i.e.
#' alpha * log(beta) - lgamma(alpha) + (alpha-1) * log(x) - beta x
#'@export
rgamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  f <- get_max_coupling(function(n) rgamma(n, alpha1, beta1),
                        function(x) dgamma(x, alpha1, beta1, log = TRUE),
                        function(n) rgamma(n, alpha2, beta2),
                        function(x) dgamma(x, alpha2, beta2, log = TRUE))
  return(f())
}
