# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = 6)

# check Gamma sampler
alpha <- 3.1
beta <- 2
x <- rgamma(100000, alpha, beta)
summary(x)
quantile(x, probs = c(0.99))
hist(x[x < 5], nclass = 200, prob = TRUE)
curve(dgamma(x, alpha, rate = beta), col = "red", add = TRUE)

### maximal coupling of Gamma
alpha1 <- 3.1
beta1 <- 2.03
alpha2 <- 7.4
beta2 <- 1.7
# sample
rgamma_coupled(alpha1, alpha2, beta1, beta2)

xy <- foreach(i = 1:50000) %dorng% {
  rgamma_coupled(alpha1, alpha2, beta1, beta2)
}

x1 <- sapply(xy, function(x) x[1])
x2 <- sapply(xy, function(x) x[2])
quantile(x1, probs = c(0.99))
hist(x1[x1 < 5], prob = TRUE, nclass = 100)
curve(dgamma(x, alpha1, rate = beta1), col = "red", add = TRUE)
quantile(x2, probs = c(0.99))
hist(sapply(xy, function(x) x[2]), prob = TRUE, nclass = 100)
curve(dgamma(x, alpha2, rate = beta2), col = "red", add = TRUE)

# same cost in both directions
microbenchmark::microbenchmark(rgamma_coupled(alpha1, alpha2, beta1, beta2),
                               rgamma_coupled(alpha2, alpha1, beta2, beta1), times = 10000)

