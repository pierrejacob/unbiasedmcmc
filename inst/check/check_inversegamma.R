# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())

# check inverse Gamma sampler
alpha <- 3.1
beta <- 2
x <- rigamma(100000, alpha, beta)
summary(x)
quantile(x, probs = c(0.99))
hist(x[x < 5], nclass = 200, prob = TRUE)
curve(exp(digamma(x, alpha, beta)), col = "red", add = TRUE)

### maximal coupling of inverse Gamma
alpha1 <- 3.1
beta1 <- 2.03
alpha2 <- 7.4
beta2 <- 1.7
# sample
rigamma_coupled(alpha1, alpha2, beta1, beta2)
# cig <- get_max_coupling(function(n) rigamma(n, alpha1, beta1),
#                  function(x) digamma(x, alpha1, beta1),
#                  function(n) rigamma(n, alpha2, beta2),
#                  function(x) digamma(x, alpha2, beta2))
xy <- foreach(i = 1:10000) %dorng% {
  rigamma_coupled(alpha1, alpha2, beta1, beta2)
}

x1 <- sapply(xy, function(x) x[1])
x2 <- sapply(xy, function(x) x[2])
quantile(x1, probs = c(0.99))
hist(x1[x1 < 5], prob = TRUE, nclass = 100)
curve(exp(digamma(x, alpha1, beta1)), col = "red", add = TRUE)
quantile(x2, probs = c(0.99))
hist(sapply(xy, function(x) x[2]), prob = TRUE, nclass = 100)
curve(exp(digamma(x, alpha2, beta2)), col = "red", add = TRUE)

# same cost in both directions
# microbenchmark::microbenchmark(rigamma_coupled(alpha1, alpha2, beta1, beta2),
#                                 cig(), times = 10000)
#
