#'@rdname fast_rmvnorm
#'@title fast_rmvnorm
#'@description Generate multivariate Normal draws. The function does not check
#' the arguments, use at your own risk.
#'@param n An integer >= 1 specifying the desired number of draws
#'@param mean A vector of size d specifying the mean vector of the multivariate Normal
#'@param covariance A matrix of size d x d specifying the covariance matrix of the multivariate Normal
#'@return A matrix of size n x d containing n d-dimensional multivariate Normal draws (one per row)
#'@examples
#'  fast_rmvnorm(2, rep(0, 5), diag(1, 5, 5))
#'@export
fast_rmvnorm <- function(n, mean, covariance){
  return(fast_rmvnorm_(n, mean, covariance))
}

#'@rdname fast_rmvnorm_chol
#'@title fast_rmvnorm_chol
#'@description Generate multivariate Normal draws. The function does not check
#' the arguments, use at your own risk.
#'@param n An integer >= 1 specifying the desired number of draws
#'@param mean A vector of size d specifying the mean vector of the multivariate Normal
#'@param chol A matrix of size d x d specifying the upper triangular Cholesky factor
#'of the covariance matrix of the multivariate Normal target,
#'for instance obtained using the \code{\link[base]{chol}}
#'function of R.
#'@return A matrix of size n x d containing n d-dimensional multivariate Normal draws (one per row)
#'@examples
#' Sigma <- diag(1, 5, 5)
#' Sigma[1,2] <- Sigma[2,1] <- 0.3
#' fast_rmvnorm_chol(2, rep(0, 5), chol(Sigma))
#'@export
fast_rmvnorm_chol <- function(nparticles, mean, chol){
  return(fast_rmvnorm_cholesky_(nparticles, mean, chol))
}

#'@rdname fast_dmvnorm
#'@title fast_dmvnorm
#'@description Compute multivariate Normal density (log-value) evaluated at each row of a given matrix. The function does not check
#' the arguments, use at your own risk.
#'@param x A matrix of size n times d
#'@param mean A vector of size d specifying the mean vector of the multivariate Normal
#'@param covariance A matrix of size d x d specifying the covariance matrix of the multivariate Normal
#'@return A vector of n evaluations of the multivariate Normal log-pdf, one for each row of \code{x}
#'@examples
#'x <- fast_rmvnorm(2, rep(0, 5), diag(1,5,5))
#'fast_dmvnorm(x, rep(0, 5), diag(1,5,5))
#'@export
fast_dmvnorm <- function(x, mean, covariance){
  return(fast_dmvnorm_(x, mean, covariance))
}

#'@rdname fast_dmvnorm_chol_inverse
#'@title fast_dmvnorm_chol_inverse
#'@description Compute multivariate Normal density (log-value) evaluated at each row of a given matrix. The function does not check
#' the arguments, use at your own risk.
#'@param x A matrix of size n times d
#'@param mean A vector of size d specifying the mean vector of the multivariate Normal
#'@param chol_inverse A matrix of size d x d specifying the inverse of the upper-triangular Cholesky
#'factor of the covariance matrix of the multivariate Normal,
#'for instance obtained using \code{solve(chol(Sigma))}
#'@return A vector of n evaluations of the multivariate Normal log-pdf, one for each row of \code{x}
#'@examples
#' Sigma <- diag(1, 5, 5)
#' Sigma[1,2] <- Sigma[2,1] <- 0.3
#' Sigma_chol <- chol(Sigma)
#' x <- fast_rmvnorm_chol(2, rep(0, 5), Sigma_chol)
#' fast_dmvnorm_chol_inverse(x, rep(0, 5), solve(Sigma_chol))
#'@export
fast_dmvnorm_chol_inverse <- function(x, mean, chol_inverse){
  return(fast_dmvnorm_cholesky_inverse_(x, mean, chol_inverse))
}



