# This script tests the functions implementing different couplings
# of Normal distributions.

# load packages
library(unbiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)


### maximal coupling of univariate Normals
mu1 <- 0.2
mu2 <- -0.8
sigma1 <- 0.4
sigma2 <- 1.7
# sample
rnorm_max_coupling(mu1, mu2, sigma1, sigma2)
# sample many times
xy <- foreach(i = 1:1000) %dorng% {
  rnorm_max_coupling(mu1, mu2, sigma1, sigma2)
}
par(mfrow = c(2,1))
hist(sapply(xy, function(x) x$xy[1]), prob = TRUE, nclass = 100)
curve(dnorm(x, mu1, sigma1), add = TRUE, col = "red")

hist(sapply(xy, function(x) x$xy[2]), prob = TRUE, nclass = 100)
curve(dnorm(x, mu2, sigma2), add = TRUE, col = "red")

mean(sapply(xy, function(x) x$identical))
PP <- integrate(f = function(x) pmin(dnorm(x, mu1, sigma1),dnorm(x, mu2, sigma2)), -20, 20)$value
PP


### maximal coupling of bivariate Gaussians
p <- 3
mu1 <- rep(0, p)
mu2 <- rep(1, p)
Sigma1 <- diag(0.4, p, p)
Sigma1[1,2] <- Sigma1[2,1] <- 0.2
Sigma2 <- diag(1.4, p, p)
Sigma2[1,2] <- Sigma2[2,1] <- -0.5
# cov(mvtnorm::rmvnorm(1e5, mu2, Sigma2))
# cov(fast_rmvnorm(1e5, mu2, Sigma2))
# colMeans(fast_rmvnorm(1e5, mu2, Sigma2))

# in "pure R"
rmvnorm_max2 <- function(mu1, mu2, Sigma1, Sigma2){
  Sigma1_chol <- chol(Sigma1)
  Sigma2_chol <- chol(Sigma2)
  Sigma1_chol_inv <- solve(Sigma1_chol)
  Sigma2_chol_inv <- solve(Sigma2_chol)
  x <- fast_rmvnorm_chol(1, mu1, Sigma1_chol)
  if (fast_dmvnorm_chol_inverse(x, mu1, Sigma1_chol_inv) + log(runif(1)) <
      fast_dmvnorm_chol_inverse(x, mu2, Sigma2_chol_inv)){
    return(list(xy = cbind(t(x),t(x)), identical = TRUE, nattempts = 0))
  } else {
    reject <- TRUE
    y <- NA
    nattempts <- 0
    while (reject){
      nattempts <- nattempts + 1
      y <- fast_rmvnorm_chol(1, mu2, Sigma2_chol)
      reject <- (fast_dmvnorm_chol_inverse(y, mu2, Sigma2_chol_inv) + log(runif(1)) <
                   fast_dmvnorm_chol_inverse(y, mu1, Sigma1_chol_inv))
    }
    return(list(xy = cbind(t(x),t(y)), identical = FALSE, nattempts = nattempts))
  }
}


# sample
nrep <- 1e4

xmax_R <- foreach (irep = 1:nrep) %dorng% {
  rmvnorm_max2(mu1, mu2, Sigma1, Sigma2)
}

xmax <- foreach (irep = 1:nrep) %dorng% {
  rmvnorm_max(mu1, mu2, Sigma1, Sigma2)
}


## check distribution of 1st sample
cat(rowMeans(sapply(xmax_R, function(x) x$xy[,1])), " versus ", mu1, "\n")
cat(rowMeans(sapply(xmax, function(x) x$xy[,1])), " versus ", mu1, "\n")

cat(cov(t(sapply(xmax_R, function(x) x$xy[,1]))), " versus\n ", Sigma1, "\n")
cat(cov(t(sapply(xmax, function(x) x$xy[,1]))), " versus\n ", Sigma1, "\n")

## check distribution of 2nd sample
cat(rowMeans(sapply(xmax_R, function(x) x$xy[,2])), " versus ", mu2, "\n")
cat(rowMeans(sapply(xmax, function(x) x$xy[,2])), " versus ", mu2, "\n")

cat(cov(t(sapply(xmax_R, function(x) x$xy[,2]))), " versus\n ", Sigma2, "\n")
cat(cov(t(sapply(xmax, function(x) x$xy[,2]))), " versus\n ", Sigma2, "\n")

microbenchmark::microbenchmark(
  a = rmvnorm_max2(mu1, mu2, Sigma1, Sigma2),
  b = rmvnorm_max(mu1, mu2, Sigma1, Sigma2),
  times = 100
)

# estimator of 1-TV:
mean(sapply(xmax, function(v) v$identical))
mean(sapply(xmax_R, function(v) v$identical))

### try the function based on Cholesky decompositions ----
Sigma1_chol <- chol(Sigma1)
Sigma1_chol_inv <- solve(chol(Sigma1))
Sigma2_chol <- chol(Sigma2)
Sigma2_chol_inv <- solve(chol(Sigma2))
xmax_chol <- foreach(i = 1:1e4) %dorng% {
  rmvnorm_max_chol(mu1, mu2, Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
}
# estimated means
cat(rowMeans(sapply(xmax_chol, function(x) x$xy[,1])), " versus ", mu1, "\n")
cat(rowMeans(sapply(xmax_chol, function(x) x$xy[,2])), " versus ", mu2, "\n")
# estimated covariances
cat(cov(t(sapply(xmax_chol, function(x) x$xy[,1]))), " versus\n", Sigma1, "\n")
cat(cov(t(sapply(xmax_chol, function(x) x$xy[,2]))), " versus\n", Sigma2, "\n")

mean(sapply(xmax_chol, function(v) v$identical))




### reflection-maximal coupling between two Normals with same covariance matrix
xreflect <- foreach(i = 1:1e4) %dorng% {
  rmvnorm_reflectionmax(mu1, mu2, Sigma1_chol, Sigma1_chol_inv)
}
# estimated means
cat(rowMeans(sapply(xreflect, function(x) x$xy[,1])), " versus ", mu1, "\n")
cat(rowMeans(sapply(xreflect, function(x) x$xy[,2])), " versus ", mu2, "\n")
# estimated covariances
cat(cov(t(sapply(xreflect, function(x) x$xy[,1]))), " versus ", Sigma1, "\n")
cat(cov(t(sapply(xreflect, function(x) x$xy[,2]))), " versus ", Sigma1, "\n")

mean(sapply(xreflect, function(v) v$identical))

## The function using Cholesky brings speedups mostly in large dimensons
p <- 500
mu1 <- rep(0, p)
mu2 <- rep(1, p)
Sigma1 <- diag(0.4, p, p)
Sigma1[1,2] <- Sigma1[2,1] <- 0.2
Sigma2 <- diag(1.4, p, p)
Sigma2[1,2] <- Sigma2[2,1] <- -0.5
Sigma1_chol <- chol(Sigma1)
Sigma1_chol_inv <- solve(chol(Sigma1))
Sigma2_chol <- chol(Sigma2)
Sigma2_chol_inv <- solve(chol(Sigma2))
microbenchmark::microbenchmark(
  a = rmvnorm_max(mu1, mu2, Sigma1, Sigma2),
  b = rmvnorm_max_chol(mu1, mu2, Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv),
  times = 100
)

