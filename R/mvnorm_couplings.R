#'@rdname rmvnorm_max
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means and covariance matrices. See \code{\link{rmvnorm_max_chol}}
#'for a version using Cholesky factors and their inverses.
#'See \code{\link{rmvnorm_reflectionmax}}
#'for a reflection-maximal coupling, in the case Sigma1=Sigma2.
#'@param mu1 First mean
#'@param mu2 First mean
#'@param Sigma1 First covariance matrix
#'@param Sigma2 Second covariance matrix
#'@return A list containing 'xy', a matrix with 2 columns (one for each draw),
#' and a boolean indicator 'identical' indicating whether the two draws
#' are identical.
#'@examples
#' p <- 3
#' mu1 <- rep(0, p)
#' mu2 <- rep(1, p)
#' Sigma1 <- diag(0.4, p, p)
#' Sigma1[1,2] <- Sigma1[2,1] <- 0.2
#' Sigma2 <- diag(1.4, p, p)
#' Sigma2[1,2] <- Sigma2[2,1] <- -0.5
#' rmvnorm_max(mu1, mu2, Sigma1, Sigma2)
#'@export
rmvnorm_max <- function(mu1, mu2, Sigma1, Sigma2){
  return(rmvnorm_max_coupling_(mu1, mu2, Sigma1, Sigma2))
}

#'@rdname rmvnorm_max_chol
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means, the Cholesky factors of their covariance matrices,
#'and the inverse of the Cholesky factors of the covariance matrices.
#'@param mu1 First mean
#'@param mu2 First mean
#'@param Cholesky1 First Cholesky factor, e.g. obtained with \code{\link[base]{chol}}
#'@param Cholesky2 Second Cholesky factor
#'@param Cholesky_inverse1 First inverse of Cholesky factor, e.g. obtained with \code{solve(chol(Sigma))}
#'@param Cholesky_inverse2 Second inverse of Cholesky factor
#'@return A list containing 'xy', a matrix with 2 columns (one for each draw),
#' and a boolean indicator 'identical' indicating whether the two draws
#' are identical.
#'@examples
#' p <- 3
#' mu1 <- rep(0, p)
#' mu2 <- rep(1, p)
#' Sigma1 <- diag(0.4, p, p)
#' Sigma1[1,2] <- Sigma1[2,1] <- 0.2
#' Sigma2 <- diag(1.4, p, p)
#' Sigma2[1,2] <- Sigma2[2,1] <- -0.5
#' Sigma1_chol <- chol(Sigma1)
#' Sigma2_chol <- chol(Sigma2)
#' Sigma1_chol_inv <- solve(Sigma1_chol)
#' Sigma2_chol_inv <- solve(Sigma2_chol)
#' rmvnorm_max_chol(mu1, mu2, Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
#'@export
rmvnorm_max_chol <- function(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2){
  # we need Cholesky <- chol(Sigma), not the transpose
  valid_Cholesky1 <- all(Cholesky1[lower.tri(Cholesky1)]==0)
  valid_Cholesky2 <- all(Cholesky2[lower.tri(Cholesky2)]==0)
  stopifnot(valid_Cholesky1, valid_Cholesky2)
  return(rmvnorm_max_coupling_cholesky(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2))
}


#'@rdname rmvnorm_reflectionmax
#'@title Reflection-Maximal coupling of two multivariate Normal distributions
#'@description Sample from reflection-maximal coupling of two multivariate Normal distributions,
#'specified through their means, with the same covariance matrix, specified
#'through its Cholesky factor and inverse of Cholesky factor.
#'
#'The idea is that a multivariate Normal is drawn around the first mean (mu1),
#'and then reflected with respect to a hyperplane orthogonal to the direction between mu1 and mu2.
#'
#'For univariate Normal distribution, see \code{\link{rnorm_reflectionmax}}.
#'
#'@param mu1 First mean
#'@param mu2 First mean
#'@param Cholesky Cholesky factor, e.g. obtained with \code{\link[base]{chol}}
#'@param Cholesky_inverse Inverse of Cholesky factor, e.g. obtained with \code{solve(chol(Sigma))}
#'@return A list containing 'xy', a matrix with 2 columns (one for each draw),
#' and a boolean indicator 'identical' indicating whether the two draws
#' are identical.
#'@examples
#' p <- 3
#' mu1 <- rep(0, p)
#' mu2 <- rep(1, p)
#' Sigma <- diag(0.4, p, p)
#' Sigma[1,2] <- Sigma[2,1] <- 0.2
#' Sigma_chol <- chol(Sigma)
#' Sigma_chol_inv <- solve(Sigma_chol)
#' rmvnorm_reflectionmax(mu1, mu2, Sigma_chol, Sigma_chol_inv)
#'@export
rmvnorm_reflectionmax <- function(mu1, mu2, Cholesky, Cholesky_inverse){
  d_ <- dim(Cholesky)
  if (is.null(d_) || (d_[1] == 1)){
    stop("function 'rmvnorm_reflectionmax' is meant for multivariate Normals; for univariate Normals, use 'rnorm_reflectionmax'")
  }
  return(rmvnorm_reflection_max_coupling_(mu1, mu2, Cholesky, Cholesky_inverse))
}


