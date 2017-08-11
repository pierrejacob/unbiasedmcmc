# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())

### maximal coupling of bivariate Gaussians
p <- 20
mu1 <- rep(0, p)
mu2 <- rep(1, p)
Sigma1 <- diag(0.4, p, p)
Sigma1[1,2] <- Sigma1[2,1] <- 0.2
Sigma2 <- diag(1.4, p, p)
Sigma2[1,2] <- Sigma2[2,1] <- -0.5
# sample

x <- fast_rmvnorm(2, mu2, Sigma2)
library(mvtnorm)
# four equivalent calculations
fast_dmvnorm(x, mu1, Sigma1)
dmvnorm(x, mu1, Sigma1, log = T)
fast_dmvnorm_chol_inverse(x, mu1, t(chol(solve(Sigma1))))
fast_dmvnorm_chol_inverse(x, mu1, solve(chol(Sigma1)))

library(microbenchmark)
microbenchmark(
  fast_dmvnorm(x, mu1, Sigma1),
  dmvnorm(x, mu1, Sigma1, log = T),
  fast_dmvnorm_chol_inverse(x, mu1, t(chol(solve(Sigma1)))),
  fast_dmvnorm_chol_inverse(x, mu1, solve(chol(Sigma1)))
)
# but if we precompute the cholesky's
chol_inv1 <- t(chol(solve(Sigma1)))
chol_inv2 <- solve(chol(Sigma1))
microbenchmark(
  fast_dmvnorm(x, mu1, Sigma1),
  dmvnorm(x, mu1, Sigma1, log = T),
  fast_dmvnorm_chol_inverse(x, mu1, chol_inv1), times = 1e3
)

# now check generator
n <- 1e4
x_1 <- fast_rmvnorm(n, mu1, Sigma1)
x_2 <- rmvnorm(n, mu1, Sigma1)
cov(x_1)[1:4,1:4]
cov(x_2)[1:4,1:4]
Sigma1[1:4,1:4]

# and using cholesky
x_3 <- fast_rmvnorm_chol(n, mu1, chol(Sigma1))
cov(x_3)[1:4,1:4]
Sigma1[1:4,1:4]
