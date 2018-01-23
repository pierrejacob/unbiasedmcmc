# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())


### maximal coupling of bivariate Gaussians
mu1 <- 0.2
mu2 <- -0.8
sigma1 <- 0.4
sigma2 <- 1.7
# sample
rnorm_max_coupling(mu1, mu2, sigma1, sigma2)

xy <- foreach(i = 1:100000) %dorng% {
  rnorm_max_coupling(mu1, mu2, sigma1, sigma2)
}
head(xy)
hist(sapply(xy, function(x) x[1]), prob = TRUE, nclass = 100)
curve(dnorm(x, mu1, sigma1), add = TRUE, col = "red")

hist(sapply(xy, function(x) x[2]), prob = TRUE, nclass = 100)
curve(dnorm(x, mu2, sigma2), add = TRUE, col = "red")

mean(sapply(xy, function(x) x[2] == x[1]))

# redefine function to count the cost
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  x <- rnorm(1, mu1, sigma1)
  if (dnorm(x, mu1, sigma1, log = TRUE) + log(runif(1)) < dnorm(x, mu2, sigma2, log = TRUE)){
    return(c(x,x,1))
  } else {
    reject <- TRUE
    y <- NA
    cost <- 1
    while (reject){
      y <- rnorm(1, mu2, sigma2)
      reject <- (dnorm(y, mu2, sigma2, log = TRUE) + log(runif(1)) < dnorm(y, mu1, sigma1, log = TRUE))
      cost <- cost + 1
    }
    return(c(x,y,cost))
  }
}

### maximal coupling of bivariate Gaussians
mu1 <- 0.2
mu2 <- -0.8
sigma1 <- 0.4
sigma2 <- 1.7
# sample
rnorm_max_coupling(mu1, mu2, sigma1, sigma2)

xy <- foreach(i = 1:10000) %dorng% {
  rnorm_max_coupling(mu1, mu2, sigma1, sigma2)
}

# estimated means
mean(sapply(xy, function(x) x[1]))
mean(sapply(xy, function(x) x[2]))
sd(sapply(xy, function(x) x[1]))
sd(sapply(xy, function(x) x[2]))
mean(sapply(xy, function(x) x[3]))
# note the same cost if we switch the distributions
xy <- foreach(i = 1:10000) %dorng% {
  rnorm_max_coupling(mu2, mu1, sigma2, sigma1)
}
mean(sapply(xy, function(x) x[3]))


### maximal coupling of bivariate Gaussians
p <- 2
mu1 <- rep(0, p)
mu2 <- rep(1, p)
Sigma1 <- diag(0.4, p, p)
Sigma1[1,2] <- Sigma1[2,1] <- 0.2
Sigma2 <- diag(1.4, p, p)
Sigma2[1,2] <- Sigma2[2,1] <- -0.5
# sample
x <- foreach(i = 1:10000) %dorng% {
  gaussian_max_coupling(mu1, mu2, Sigma1, Sigma2)
}
# estimated means
rowMeans(sapply(x, function(x) x[,1]))
rowMeans(sapply(x, function(x) x[,2]))
# estimated covariances
cov(t(sapply(x, function(x) x[,1])))
cov(t(sapply(x, function(x) x[,2])))
### also try the function based on Cholesky decompositions


Sigma1_chol <- chol(Sigma1)
Sigma1_chol_inv <- solve(chol(Sigma1))
Sigma2_chol <- chol(Sigma2)
Sigma2_chol_inv <- solve(chol(Sigma2))

x <- foreach(i = 1:50000) %dorng% {
  gaussian_max_coupling_cholesky_R(mu1, mu2, Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
}
# estimated means
rowMeans(sapply(x, function(x) x[,1]))
rowMeans(sapply(x, function(x) x[,2]))
# estimated covariances
cov(t(sapply(x, function(x) x[,1])))
cov(t(sapply(x, function(x) x[,2])))


### optimal transport coupling of bivariate Gaussians
p <- 3
mu1 <- rep(0, p)
mu2 <- rep(1, p)
Sigma1 <- diag(0.4, p, p)
Sigma1[1,2] <- Sigma1[2,1] <- 0.2
Sigma2 <- diag(1.4, p, p)
Sigma2[1,2] <- Sigma2[2,1] <- -0.5

S1_chol <- chol(Sigma1)
S1_chol_inv <- solve(chol(Sigma1))
S2_chol <- chol(Sigma2)
S2_chol_inv <- solve(chol(Sigma2))

#### sample
x <- gaussian_opt_transport(10000, mu1, mu2, S1_chol, S2_chol, S1_chol_inv, S2_chol_inv)

rowMeans(sapply(x, function(x) x[,1]))
rowMeans(sapply(x, function(x) x[,2]))

cov(t(sapply(x, function(x) x[,1])))
cov(t(sapply(x, function(x) x[,2])))







