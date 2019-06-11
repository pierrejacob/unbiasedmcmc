#'@rdname dinvgaussian
#'@title Log-density of inverse Gaussian
#'@description Computes log-pdf at x > 0 of inverse Gaussian with given parameters mu, lambda, given by
#'  \deqn{0.5 * log(\lambda/(2*\pi)) - 1.5 * log(x) - \lambda * (x-\mu)^2 / (2 * \mu^2 * x)}
#'@return A vector of n log-pdf values, one for each element in the first argument 'x'.
#'@export
dinvgaussian <- function(x, mu, lambda){
  return(0.5 * log(lambda/(2*pi)) - 1.5 * log(x) - lambda * (x-mu)^2 / (2 * mu^2 * x))
}

#'@rdname rinvgaussian
#'@title Sample from inverse Gaussian
#'@description Parametrized by mu, lambda, with log-density given by
#'  \deqn{0.5 * log(\lambda/(2*\pi)) - 1.5 * log(x) - \lambda * (x-\mu)^2 / (2 * \mu^2 * x)}
#'
#'  The procedure goes as follows.
#'
#'  \itemize{
#'  \item Generate nu ~ Normal(0,1).
#'  \item Define y = nu^2.
#'  \item Define x = mu + mu^2 * y / (2 * lambda) - mu / (2 * lambda) * sqrt(4 * mu * lambda * y + mu^2 * y^2).
#'  \item Generate Z ~ Uniform(0,1).
#'  \item If z <= mu / (mu + x), output x, otherwise output mu^2 / x.
#'  }
#'@return A vector of n draws, where n is the first argument.
#'@export
rinvgaussian <- function(n, mu, lambda){
  return(rinvgaussian_c(n, mu, lambda))
}

#'@rdname rinvgaussian_coupled
#'@title Sample from maximally coupled inverse Gaussian
#'@description with parameters mu1, mu2, lambda1, lambda2; see \code{\link{rinvgaussian}}.
#'@return A pair of values in a vector of size two.
#'@export
rinvgaussian_coupled <- function(mu1, mu2, lambda1, lambda2){
  rinvgaussian_coupled_c(mu1, mu2, lambda1, lambda2)
}
