
# from gaussian_max_coupling ----------------------------------------------


#'@export
gaussian_max_coupling_cholesky_R <- function(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2){
  # we need Cholesky <- chol(Sigma), not the transpose
  valid_Cholesky1 <- all(Cholesky1[lower.tri(Cholesky1)]==0)
  valid_Cholesky2 <- all(Cholesky2[lower.tri(Cholesky2)]==0)
  stopifnot(valid_Cholesky1, valid_Cholesky2)

  return(gaussian_max_coupling_cholesky(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2))
}

#'@export
gaussian_max_coupling <- function(mu1, mu2, Sigma1, Sigma2){
  return(gaussian_max_couplingC(mu1, mu2, Sigma1, Sigma2))
}

# from gaussian_opt_transport ---------------------------------------------

#'@export
gaussian_opt_transport <- function(nsamples, m1, m2,
                                   S1_chol,
                                   S2_chol,
                                   S1_chol_inv,
                                   S2_chol_inv){
  # the cholesky factors should be what comes out of chol(.), not the transpose
  valid_S1_chol <- all(S1_chol[lower.tri(S1_chol)]==0)
  valid_S2_chol <- all(S2_chol[lower.tri(S2_chol)]==0)
  stopifnot(valid_S1_chol, valid_S2_chol)

  #S2 <- crossprod(S2_chol)
  #transport_matrix <- S1_chol_inv %*% chol(S1_chol %*% S2 %*% t(S1_chol)) %*% t(S1_chol_inv)

  # this is equivalent to but a bit faster than the above
  transport_matrix <- S1_chol_inv %*% chol(tcrossprod(tcrossprod(S1_chol,S2_chol))) %*% t(S1_chol_inv)

  x <- vector(nsamples,mode='list')
  m1 <- t(as.matrix(m1))
  m2 <- t(as.matrix(m2))

  # m1 and m2 need to come in as vectors or p x 1 matrices so that they
  # are 1 x p matrices at this point in the code
  valid_m1 <- all(dim(m1) == c(1,nrow(S1_chol)))
  valid_m2 <- all(dim(m2) == c(1,nrow(S2_chol)))
  stopifnot(valid_S1_chol, valid_S2_chol)

  for (i in 1:nsamples){
    x1 <- fast_rmvnorm_chol(1,m1,S1_chol)
    x2 <- (sweep(x1, 2, m1, FUN = "-")) %*% transport_matrix
    x2 <- sweep(x2, 2, m2, FUN = "+")
    x[[i]] <- cbind.data.frame(x1=drop(x1),x2=drop(x2))
  }
  return(x)
}

