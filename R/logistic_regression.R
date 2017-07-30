
# This file consists of functions related estimating logistic regression models using the original
# PGG algorithm and our unbiased estimator

# from util_logisticregression --------------------------------------------
#'@export
xbeta <- function(X, beta){
  return(xbeta_(X, beta))
}

#'@export
sigma <- function(X, w){
  return(sigma_(X, w))
}

#'@export
sigma_function <- function(omega, X, invB){
  return(solve(sigma(X, omega) + invB))
}

#'@export
m_function <- function(omega, Sigma, X, Y, invBtimesb){
  return(Sigma %*% (t(X) %*% matrix(Y - rep(0.5, length(Y)), ncol = 1) + invBtimesb))
}

# from logistic_precomputation --------------------------------------------

#'@rdname logistic_precomputation
#'@title Precomputation to prepare for the Polya-Gamma sampler
#'@description This function takes the canonical elements defining the logistic regression
#' problem (the vector of outcome Y, covariate matrix X, the prior mean b and the prior variance B),
#' and precomputes some quantities repeatedly used in the Polya-Gamma sampler and variants of it.
#' The precomputed quantities are returned in a list, which is then meant to be passed to the samplers.
#'@export
logistic_precomputation <- function(Y, X, b, B){
  invB <- solve(B)
  invBtimesb <- invB %*% b
  Ykappa <- matrix(Y - rep(0.5, length(Y)), ncol=1)
  XTkappa <- t(X) %*% Ykappa
  KTkappaplusinvBtimesb <- XTkappa + invBtimesb
  return(list(n=nrow(X), p=ncol(X), X=X, Y=Y, b=b, B=B,
              invB=invB, invBtimesb=invBtimesb, KTkappaplusinvBtimesb=KTkappaplusinvBtimesb))
}


# from m_and_sigma --------------------------------------------------------

# The following function computes m(omega) and Sigma(omega)... (or what we really need instead)
# it returns m (= m(omega)), Sigma_inverse = Sigma(omega)^{-1},
# as well as Cholesky_inverse and Cholesky that are such that
# Cholesky_inverse is the lower triangular matrix L, in the decomposition Sigma^{-1} = L L^T
# whereas Cholesky is the lower triangular matrix Ltilde in the decomposition Sigma = Ltilde^T Ltilde
#'@export
m_and_sigma <- function(omega, X, invB, KTkappaplusinvBtimesb){
  return(m_sigma_function_(omega, X, invB, KTkappaplusinvBtimesb))
}


# from pg_gibbs -----------------------------------------------------------
#'@rdname pg_gibbs
#'@title Polya-Gamma Gibbs sampler
#'@description This implements the sampler proposed in
#' Nicholas G Polson, James G Scott, and Jesse Windle. Bayesian inference for logistic models using
#' Polya–Gamma latent variables. Journal of the American statistical Association, 108(504):1339–1349, 2013.
#' The arguments are:
#'  \itemize{
#'   \item niterations: the number of desired MCMC iterations,
#'   \item logistic_setting: a list of precomputed quantities obtained via 'logistic_precomputation'.
#'   }
#'@return a matrix where each row corresponds to an iteration of the sampler, and contains
#' the regression coefficient at that iteration.
#'@export
pg_gibbs <- function(niterations, logistic_setting){
  X <- logistic_setting$X
  # Y <- logistic_setting$Y
  # invBtimesb <- logistic_setting$invBtimesb
  n <- nrow(X)
  p <- ncol(X)
  betas <- matrix(0, ncol=p, nrow=niterations)
  beta <- matrix(0, ncol=p, nrow=1)
  for (iteration in 1:niterations){
    zs <- abs(xbeta(X, beta))
    w <- rpg(n, h=1, z=zs)
    res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
    # S <- sigma_function(w, X, invB)
    # m <- m_function(w, S, X, Y, invBtimesb)
    beta <- fast_rmvnorm_chol(1, res$m, res$Cholesky)
    # beta <- fast_rmvnorm(1, m, S)
    betas[iteration,] <- beta
  }
  return(betas)
}

# from sample_w -----------------------------------------------------------
w_indep <- function(beta1,beta2,X){
  n <- nrow(X)
  w1 <- rep(0., n)
  w2 <- rep(0., n)
  z1s <- abs(xbeta(X, beta1))
  z2s <- abs(xbeta(X, beta2))

  for (i in 1:n){
    w1[[i]] <- rpg(num=1, h=1, z=z1s[i])
    w2[[i]] <- rpg(num=1, h=1, z=z2s[i])
  }
  return(list(w1=w1, w2=w2))
}

#'@export
w_rejsampler_caller <- function(beta1,beta2,X){
  w_rejsamplerC(beta1, beta2, X)
}

#'@export
w_rejsampler <- function(beta1, beta2, X){
  # we can do rejection sampling, sampling from z1 and aiming for z2,
  # provided z2 > z1. The ratio of densities is proportional to
  logratio <- function(omega, z_min, z_max){
    return(-omega * 0.5 * (z_max^2 - z_min^2))
  }

  w1 <- rep(0., n)
  w2 <- rep(0., n)
  z1s <- abs(xbeta(X, beta1))
  z2s <- abs(xbeta(X, beta2))

  for (i in 1:n){
    z_i <- c(z1s[i], z2s[i])
    z_min <- min(z_i)
    z_max <- max(z_i)

    proposal <- rpg(num=1, h=1, z=z_min)
    w_min <- proposal

    if (log(runif(1)) < logratio(proposal, z_min, z_max)){
      w_max <- proposal
    } else {
      w_max <- rpg(num=1, h=1, z=z_max)
    }

    if (which.min(z_i) == 1){
      w1[i] <- w_min
      w2[i] <- w_max
    } else {
      w2[i] <- w_min
      w1[i] <- w_max
    }
  }
  return(list(w1=w1, w2=w2))
}

w_max_coupling <- function(beta1, beta2, X){
  w1s <- rep(0., n)
  w2s <- rep(0., n)
  z1s <- abs(xbeta(X, beta1))
  z2s <- abs(xbeta(X, beta2))

  for (i in 1:n){
    z1 <- z1s[i]
    z2 <- z2s[i]

    w1 <- rpg(num=1, h=1, z=z1)
    w1s[[i]] <- w1
    u <- runif(1,0,cosh(z1/2)*exp(-0.5*z1^2*w1))
    if(u <= cosh(z2/2)*exp(-0.5*z2^2*w1)){
      w2 <- w1
    } else {
      accept <- FALSE
      while(!accept){
        w2 <- rpg(num=1, h=1, z=z2)
        u <- runif(1,0,cosh(z2/2)*exp(-0.5*z2^2*w2))
        if(u > cosh(z1/2)*exp(-0.5*z1^2*w2)){
          accept <- TRUE
        }
      }
    }
    w2s[[i]] <- w2
  }
  return(list(w1=w1s, w2=w2s))
}

#'@export
sample_w <- function(beta1, beta2, X, mode='rej_samp'){
  if (mode == 'rej_samp'){
    w1w2_mat <- w_rejsamplerC(beta1, beta2, X)
    w1w2 <- list(w1=w1w2_mat[,1], w2=w1w2_mat[,2])
  } else if (mode=='indep'){
    w1w2 <- w_indep(beta1,beta2,X)
  } else if (mode=='max'){
    w1w2 <- w_max_coupling(beta1,beta2,X)
  }
  return(w1w2)
}


# from sample_beta --------------------------------------------------------
#'@export
sample_beta <- function(w1, w2, logistic_setting, mode="max_coupling", mc_prob=0.5){
  X <- logistic_setting$X
  KTkappaplusinvBtimesb <- logistic_setting$KTkappaplusinvBtimesb
  invB <- logistic_setting$invB
  res1 <- m_and_sigma(w1, X, invB, KTkappaplusinvBtimesb)
  res2 <- m_and_sigma(w2, X, invB, KTkappaplusinvBtimesb)

  if(mode=='max_coupling'){
    x <- gaussian_max_coupling_cholesky(res1$m, res2$m, res1$Cholesky, res2$Cholesky, res1$Cholesky_inverse, res2$Cholesky_inverse)
    beta1 <- x[,1]
    beta2 <- x[,2]
  } else if(mode=='opt_transport'){
    x <- gaussian_opt_transport(1, res1$m, res2$m, res1$Cholesky, res2$Cholesky, res1$Cholesky_inverse, res2$Cholesky_inverse)
    beta1 <- x[[1]][,1]
    beta2 <- x[[1]][,2]
  } else if(mode=='both'){
    if (runif(1) < mc_prob){
      x <- gaussian_max_coupling_cholesky(res1$m, res2$m, res1$Cholesky, res2$Cholesky, res1$Cholesky_inverse, res2$Cholesky_inverse)
      beta1 <- x[,1]
      beta2 <- x[,2]
    } else {
      x <- gaussian_opt_transport(1, res1$m, res2$m, res1$Cholesky, res2$Cholesky, res1$Cholesky_inverse, res2$Cholesky_inverse)
      beta1 <- x[[1]][,1]
      beta2 <- x[[1]][,2]
    }
  } else {
    stop('invalid coupling method')
  }
  return(list(beta1=beta1, beta2=beta2))
}

