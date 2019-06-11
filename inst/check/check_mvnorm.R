# This script tests the functions implementing different
# random number generators and probability density functions for
# multivariate Normal distributions.

# load packages
library(debiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

### sampling trivariate Gaussians
p <- 3
mu1 <- rep(0, p)
mu2 <- rep(1, p)
Sigma1 <- diag(0.4, p, p)
Sigma1[1,2] <- Sigma1[2,1] <- 0.2
Sigma2 <- diag(1.4, p, p)
Sigma2[1,2] <- Sigma2[2,1] <- -0.5

# Compute Cholesky factors and inverses
Sigma1_chol <- chol(Sigma1)
Sigma2_chol <- chol(Sigma2)
Sigma1_chol_inv <- solve(Sigma1_chol)
Sigma2_chol_inv <- solve(Sigma2_chol)
## Note on Cholesky factorization:
## we have
# t(Sigma1_chol) %*% Sigma1_chol
# and Sigma1_chol is upper triangular

# sample
x <- fast_rmvnorm(2, mu2, Sigma2)

# equivalent calculations of pdf evaluations
fast_dmvnorm(x, mu1, Sigma1)
mvtnorm::dmvnorm(x, mu1, Sigma1, log = T)
fast_dmvnorm_chol_inverse(x, mu1, Sigma1_chol_inv)

# and again
x <- fast_rmvnorm(10, mu1, Sigma1)
fast_dmvnorm(x, mu2, Sigma2)
mvtnorm::dmvnorm(x, mu2, Sigma2, log = T)
fast_dmvnorm_chol_inverse(x, mu2, Sigma2_chol_inv)

# benchmark
microbenchmark::microbenchmark(
  fast_dmvnorm(x, mu1, Sigma1),
  mvtnorm::dmvnorm(x, mu1, Sigma1, log = T),
  fast_dmvnorm_chol_inverse(x, mu1, Sigma1_chol_inv))

# now check random number generator
n <- 1e4
x2_1 <- fast_rmvnorm(n, mu2, Sigma2)
x2_2 <- mvtnorm::rmvnorm(n, mu2, Sigma2)

# compare means
cat(colMeans(x2_1), "\n", colMeans(x2_2), "\n", mu2, "\n")

# compare covariances
cat(cov(x2_1)[1:2,1:2], "\n", cov(x2_2)[1:2,1:2], "\n", Sigma2[1:2,1:2], "\n")


# and using cholesky
x1_1 <- fast_rmvnorm(n, mu1, Sigma1)
x1_3 <- fast_rmvnorm_chol(n, mu1, Sigma1_chol)

# compare means
cat(colMeans(x1_1), "\n", colMeans(x1_3), "\n", mu1, "\n")

# compare covariances
cat(cov(x1_1)[1:3,1:3], "\n", cov(x1_3)[1:3,1:3], "\n", Sigma1[1:3,1:3], "\n")
