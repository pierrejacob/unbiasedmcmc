# load packages
library(debiasedmcmc)
library(ggthemes)
library(coda)
setwd("~/Dropbox/UnbiasedMCMCResults/")
# library(parallel)
#
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = 10)

data(diabetes)
X <- scale(diabetes$x)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)
lambda <- 1
pb <- get_blasso(Y, X, lambda)
p <- ncol(X)
n <- nrow(X)
XtX <- t(X) %*% X
XtY <- t(X) %*% Y
alpha1 <- (n-1)/2 + p/2
lambda2 <- lambda^2


gibbs_kernel <- function(current_state, ...){
  beta <- current_state[1:p]
  tau2 <- current_state[(p+1):(2*p)]
  sigma2 <- current_state[2*p+1]
  res_ <- debiasedmcmc:::blassoconditional(Y, X, XtY, XtX, tau2, sigma2)
  # D_tau_inv <- diag(1/tau2, p, p)
  # A <- XtX + D_tau_inv
  # A_inv <- solve(A)
  # # update beta
  # beta <- t(fast_rmvnorm(1, (A_inv %*% XtY)[,1], sigma2 * A_inv))
  # update sigma
  # norm <- sum((Y - X %*% beta)^2)
  # betaDbeta <- sum(beta^2 / tau2)
  beta <- res_$beta
  norm <- res_$norm
  betaDbeta <- res_$betaDbeta
  sigma2 <- rigamma(1, alpha1, 0.5 * (norm + betaDbeta))
  # update tau
  sqrtlambda2sigma2 <- sqrt(lambda2 * sigma2)
  for (component in 1:p){
    tau2[component] <- 1 / rinvgaussian(1, sqrtlambda2sigma2 / abs(beta[component]), lambda2)
  }
  return(c(beta, tau2, sigma2))
}


state <- pb$rinit()
for (i in 1:10){
  state <- gibbs_kernel(state)
}
state
