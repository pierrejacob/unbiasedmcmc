library(unbiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

# check inverse Gamma sampler
alpha <- 3.1
beta <- 2
x <- rinversegamma(100000, alpha, beta)
summary(x)
quantile(x, probs = c(0.99))
hist(x[x < 5], nclass = 200, prob = TRUE)
curve(exp(dinversegamma(x, alpha, beta)), col = "red", add = TRUE)

### maximal coupling of inverse Gamma
alpha1 <- 3.1
beta1 <- 2.03
alpha2 <- 7.4
beta2 <- 1.7
# sample
rinversegamma_coupled(alpha1, alpha2, beta1, beta2)
xy <- foreach(i = 1:10000) %dorng% {
  rinversegamma_coupled(alpha1, alpha2, beta1, beta2)
}

x1 <- sapply(xy, function(x) x$xy[1])
x2 <- sapply(xy, function(x) x$xy[2])
quantile(x1, probs = c(0.99))

hist(x1[x1 < 5], prob = TRUE, nclass = 100)
curve(exp(dinversegamma(x, alpha1, beta1)), col = "red", add = TRUE)

hist(x2, prob = TRUE, nclass = 100)
curve(exp(dinversegamma(x, alpha2, beta2)), col = "red", add = TRUE)
