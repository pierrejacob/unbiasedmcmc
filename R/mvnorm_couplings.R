#'@rdname rmvnorm_max
#'@title Maximal coupling of two multivariate Normal distributions
#'@description sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means and covariance matrices
#'@export
rmvnorm_max <- function(mu1, mu2, Sigma1, Sigma2){
  return(rmvnorm_max_coupling_(mu1, mu2, Sigma1, Sigma2))
}

#'@rdname rmvnorm_max_chol
#'@title Maximal coupling of two multivariate Normal distributions
#'@description sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means, the Cholesky factors of their covariance matrices,
#'and the Cholesky factors of the inverse covariance matrices (i.e. the precision matrices)
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
#'@description sample from reflection-maximal coupling of two multivariate Normal distributions,
#'specified through their means, with the same covariance matrix specified
#'through its Cholesky factor and Cholesky factor of its inverse.
#'@export
rmvnorm_reflectionmax <- function(mu1, mu2, Cholesky, Cholesky_inverse){
  return(rmvnorm_reflection_max_coupling_(mu1, mu2, Cholesky, Cholesky_inverse))
}


