#'@rdname dinvgaussian
#'@title compute log-density of inverse Gaussian
#'@description at x > 0 and with given parameters mu, lambda, given by
#'  0.5 * log(lambda/(2*pi)) - 1.5 * log(x) - lambda * (x-mu)^2 / (2 * mu^2 * x)
#'@export
dinvgaussian <- function(x, mu, lambda){
  return(0.5 * log(lambda/(2*pi)) - 1.5 * log(x) - lambda * (x-mu)^2 / (2 * mu^2 * x))
}

#'@rdname rigamma
#'@title Sample from inverse Gaussian
#'@description with given parameters mu, lambda, with log-density given by
#'  0.5 * log(lambda/(2*pi)) - 1.5 * log(x) - lambda * (x-mu)^2 / (2 * mu^2 * x).
#'  The procedure goes as follows.
#'  Generate nu ~ Normal(0,1).
#'  Define y = nu^2.
#'  Define x = mu + mu^2 * y / (2 * lambda) - mu / (2 * lambda) * sqrt(4 * mu * lambda * y + mu^2 * y^2).
#'  Generate Z ~ Uniform(0,1).
#'  If z <= mu / (mu + x), output x, otherwise output mu^2 / x.
#'@export
rinvgaussian <- function(n, mu, lambda){
  return(rinvgaussian_c(n, mu, lambda))
}

#'@rdname rinvgaussian_coupled
#'@title Sample from maximally coupled inverse Gaussian
#'@description with given parameters mu1, mu2, lambda1, lambda2
#'@export
rinvgaussian_coupled <- function(mu1, mu2, lambda1, lambda2){
  rinvgaussian_coupled_c(mu1, mu2, lambda1, lambda2)
}
