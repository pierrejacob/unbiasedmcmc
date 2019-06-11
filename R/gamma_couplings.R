#'@rdname rgamma_coupled
#'@title Sample from maximally coupled Gamma variables
#'@description Draws a pair of variables, respectively Gamma(alpha1, beta1) and Gamma(alpha2, beta2)
#'where the parametrization is that beta is the rate, i.e. the log-pdf of Gamma(alpha,beta) evaluated at x is
#' \deqn{\alpha * log(\beta) - lgamma(\alpha) + (\alpha-1) * log(x) - \beta x}
#' where \eqn{lgamma} stands for the logarithm of the Gamma function.
#'@param alpha1 First shape
#'@param alpha2 Second shape
#'@param beta1 First rate
#'@param beta2 Second rate
#'@return A list with entry 'xy' for the pair of values, and boolean 'identical' indicating whether the two values
#'are identical.
#'@export
rgamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  f <- get_max_coupling(function(n) rgamma(n, alpha1, beta1),
                        function(x) dgamma(x, alpha1, beta1, log = TRUE),
                        function(n) rgamma(n, alpha2, beta2),
                        function(x) dgamma(x, alpha2, beta2, log = TRUE))
  return(f())
}
