#'@rdname fast_rmvnorm
#'@title fast_rmvnorm
#'@description fast_rmvnorm
#'@export
fast_rmvnorm <- function(nparticles, mean, covariance){
  return(fast_rmvnorm_(nparticles, mean, covariance))
}

#'@export
fast_rmvnorm_chol <- function(nparticles, mean, chol){
  return(fast_rmvnorm_cholesky_(nparticles, mean, chol))
}

#'@rdname fast_dmvnorm
#'@title fast_dmvnorm
#'@description fast_dmvnorm
#'@export
fast_dmvnorm <- function(x, mean, covariance){
  return(fast_dmvnorm_(x, mean, covariance))
}

#'@export
fast_dmvnorm_chol_inverse <- function(x, mean, chol_inverse){
  return(fast_dmvnorm_cholesky_inverse_(x, mean, chol_inverse))
}



