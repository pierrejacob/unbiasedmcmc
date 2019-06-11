library(debiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

### maximal coupling of bivariate Gaussians with identical covariance matrices
mu1 <- c(0.2, 0.3)
mu2 <- c(0.0, 0.8)
Sigma <- diag(1, 2, 2)
Sigma[1,2] <- Sigma[2,1] <- 0.8
Sigma[2,2] <- 2
Sigma_chol <- chol(Sigma)
Sigma_inv_chol <- solve(Sigma_chol)

# sample once
rmvnorm_reflectionmax(mu1, mu2, Sigma_chol, Sigma_inv_chol)
# function output samples (column-bound) in $xy, and a boolean indicator of identity in $identical

nsamples <- 5e4
xy <- foreach(i = 1:nsamples) %dorng% {
  rmvnorm_reflectionmax(mu1, mu2, Sigma_chol, Sigma_inv_chol)
}

identicals <- sapply(xy, function(x) x$identical)
equalvalues <- sapply(xy, function(x) all(x$xy[,1] == x$xy[,2]))
all(identicals == equalvalues)

# collect samples from both distributions
sample1 <- t(sapply(xy, function(x) x$xy[,1]))
sample2 <- t(sapply(xy, function(x) x$xy[,2]))

# and check that they follow the correct distribution
colMeans(sample1)
mu1
colMeans(sample2)
mu2

cov(sample1)
cov(sample2)
Sigma

# estimate of 1-TVD
mean(sapply(xy, function(x) x$identical))

# visualize marginals
hist(sample1[,1], prob = TRUE, nclass = 100)
curve(dnorm(x, mu1[1], sqrt(Sigma[1,1])), add = TRUE, col = "red")

hist(sample1[,2], prob = TRUE, nclass = 100)
curve(dnorm(x, mu1[2], sqrt(Sigma[2,2])), add = TRUE, col = "red")

hist(sample2[,1], prob = TRUE, nclass = 100)
curve(dnorm(x, mu2[1], sqrt(Sigma[1,1])), add = TRUE, col = "red")

hist(sample2[,2], prob = TRUE, nclass = 100)
curve(dnorm(x, mu2[2], sqrt(Sigma[2,2])), add = TRUE, col = "red")



mu1 <- 1
mu2 <- 2
Sigma_proposal <- diag(2, 1, 1)
Sigma_chol <- chol(Sigma_proposal)
Sigma_chol_inv <- solve(Sigma_chol)
nsamples <- 5e4
xy <- foreach(i = 1:nsamples) %dorng% {
  rmvnorm_reflectionmax(mu1, mu2, Sigma_chol, Sigma_inv_chol)
}
hist(sapply(xy, function(x) x$xy[,1]), prob = TRUE, nclass = 100)
curve(dnorm(x, mu1, Sigma_chol[1]), add = T)

hist(sapply(xy, function(x) x$xy[,2]), prob = TRUE, nclass = 100)
curve(dnorm(x, mu2, Sigma_chol[1]), add = T)
