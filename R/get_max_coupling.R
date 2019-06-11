#'@rdname get_max_coupling
#'@title Sample from maximal coupling of two distributions p and q
#'@description Takes two univariate continuous distributions (specified by random number generator and log-pdf function),
#' and returns a function to sample from a maximal coupling of these two distributions.
#'@param rp A function taking n as an argument and returning n samples from the distribution p
#'@param dp A function taking x as an argument and returning log-pdf of p evaluated at x
#'@param rq A function taking n as an argument and returning n samples from the distribution q
#'@param dq A function taking x as an argument and returning log-pdf of q evaluated at x
#'@return Returns a list with
#'
#' \itemize{
#' \item "xy": the pair of samples \eqn{(x,y)}
#'
#' \item "identical": TRUE if \eqn{x = y}, FALSE otherwise
#' }
#'@examples
#' mu1 <- 0; mu2 <- 1; sigma1 <- 0.5; sigma2 <- 1.2
#'  f <- get_max_coupling(function(n) rnorm(n, mu1, sigma1),
#'  function(x) dnorm(x, mu1, sigma1, log = TRUE),
#'  function(n) rnorm(n, mu2, sigma2),
#'  function(x) dnorm(x, mu2, sigma2, log = TRUE))
#' f()
#'@export
get_max_coupling <- function(rp, dp, rq, dq){
  function(){
    x <- rp(1)
    if (dp(x) + log(runif(1)) < dq(x)){
      return(list(xy = c(x,x), identical = TRUE))
    } else {
      reject <- TRUE
      y <- NA
      while (reject){
        y <- rq(1)
        reject <- (dq(y) + log(runif(1)) < dp(y))
      }
      return(list(xy = c(x,y), identical = FALSE))
    }
  }
}
