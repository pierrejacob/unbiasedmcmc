#'@rdname rnorm_max_coupling
#'@title Maximal coupling of two univariate Normal distributions
#'@description Sample from maximal coupling of two univariate Normal distributions,
#'specified through their means and standard deviations. See \code{\link{rmvnorm_max}} for a multivariate version.
#'@param mu1 First mean
#'@param mu2 Second mean
#'@param sigma1 First mean
#'@param sigma2 Second mean
#'@return Returns a list with
#'
#' \itemize{
#' \item "xy": the pair of samples \eqn{(x,y)}
#'
#' \item "identical": TRUE if \eqn{x = y}, FALSE otherwise
#' }
#'@export
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  f <- get_max_coupling(function(n) rnorm(n, mu1, sigma1),
                        function(x) dnorm(x, mu1, sigma1, log = TRUE),
                        function(n) rnorm(n, mu2, sigma2),
                        function(x) dnorm(x, mu2, sigma2, log = TRUE))
  return(f())
}
