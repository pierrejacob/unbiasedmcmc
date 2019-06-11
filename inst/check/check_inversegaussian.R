library(debiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

# check inverse Gamma sampler
mu <- 3.1
lambda <- 2
x <- rinvgaussian(100000, mu, lambda)
summary(x)
quantile(x, probs = c(0.99))
hist(x, nclass = 200, prob = TRUE)
curve(exp(dinvgaussian(x, mu, lambda)), col = "red", add = TRUE)

### maximal coupling of inverse Gamma
mu1 <- 3.1
lambda1 <- 2
mu2 <- 1.7
lambda2 <- 1.3
# sample
rinvgaussian_coupled(mu1, mu2, lambda1, lambda2)

xy <- foreach(i = 1:10000) %dorng% {
  rinvgaussian_coupled(mu1, mu2, lambda1, lambda2)
}

x1 <- sapply(xy, function(x) x[1])
x2 <- sapply(xy, function(x) x[2])
quantile(x1, probs = c(0.99))
hist(x1, prob = TRUE, nclass = 100)
curve(exp(dinvgaussian(x, mu1, lambda1)), col = "red", add = TRUE)
quantile(x2, probs = c(0.99))

hist(x2, prob = TRUE, nclass = 100)
curve(exp(dinvgaussian(x, mu2, lambda2)), col = "red", add = TRUE)
