# load packages
library(debiasedmcmc)
setwd("~/Dropbox/UnbiasedMCMCResults/")
rm(list = ls())
set.seed(1)

## regression setting
n <- 100
p <- 500
SNR <- 2
s_star <- 10
s0 <- 100
sigma0 <- 1
# beta_star <- SNR * sqrt(sigma0^2 * log(p) / n) * c(2,-3,2,2,-3,3,-2,3,-2,3, rep(0, p-10))
beta_star <- rnorm(p)
# independent design
X <- matrix(rnorm(n * p), nrow = n, ncol = p) # fast_rmvnorm_chol(n, rep(0, p), diag(1, p, p))
X <- scale(X)
Y <- X %*% matrix(beta_star, ncol = 1) + rnorm(n, 0, sigma0)
Y <- scale(Y)
##

Phi <- X / sigma0
alpha <- Y / sigma0
D <- sigma0^2 * diag(1, p,p)

Sigma <- solve(t(Phi) %*% Phi + solve(D))
mu <- Sigma %*% t(Phi) %*% alpha

nsamples <- 10000
xsamples <- fast_rmvnorm(nsamples, mu[,1], Sigma)
cov(xsamples)
Sigma
colMeans(xsamples)
t(mu)

# Bhattacharya et al
PhiDPhi <- Phi %*% D %*% t(Phi)
DPhi <- D %*% t(Phi)
In <- diag(1, n, n)
bh <- function(){
  u <- fast_rmvnorm(1, rep(0, p), D)
  delta <- fast_rmvnorm(1, rep(0, n), In)
  nu <- Phi %*% t(u) + t(delta)
  w <- solve(a = PhiDPhi + In, b = alpha - nu)
  theta <- u + t(DPhi %*% w)
  return(theta)
}

registerDoParallel(cores = detectCores())
xsamples2 <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
  bh()
}

cov(xsamples2)
Sigma
colMeans(xsamples2)
t(mu)
