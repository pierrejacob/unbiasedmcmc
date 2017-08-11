# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = 6)

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
hist(sapply(xy, function(x) x[2]), prob = TRUE, nclass = 100)
curve(exp(dinvgaussian(x, mu2, lambda2)), col = "red", add = TRUE)

# same cost in both directions
microbenchmark::microbenchmark(rinvgaussian_coupled(mu1, mu2, lambda1, lambda2),
                               rinvgaussian_coupled(mu2, mu1, lambda2, lambda1), times = 10000)


xy <- foreach(i = 1:10000) %dorng% {
  rinvgaussian_coupled_2(mu1, mu2, lambda1, lambda2)
}

x1 <- sapply(xy, function(x) x[1])
x2 <- sapply(xy, function(x) x[2])
quantile(x1, probs = c(0.99))
hist(x1, prob = TRUE, nclass = 100)
curve(exp(dinvgaussian(x, mu1, lambda1)), col = "red", add = TRUE)
quantile(x2, probs = c(0.99))
hist(sapply(xy, function(x) x[2]), prob = TRUE, nclass = 100)
curve(exp(dinvgaussian(x, mu2, lambda2)), col = "red", add = TRUE)

# same cost in both directions
microbenchmark::microbenchmark(rinvgaussian_coupled_2(mu1, mu2, lambda1, lambda2),
                               rinvgaussian_coupled_2(mu2, mu1, lambda2, lambda1), times = 10000)



# ## Faster implementation in C++
# cppFunction("NumericVector rinvgaussian_c(int n, double mu, double lambda){
#   NumericVector results(n);
#   NumericVector nu = rnorm(n);
#   NumericVector z = runif(n);
#   NumericVector x = nu * nu;
#   x = mu + mu*mu * x / (2. * lambda) - mu / (2. * lambda) * sqrt(4. * mu * lambda * x + mu*mu * x*x);
#   for (int i = 0; i < n; i++){
#     if (z(i) <= mu / (mu + x(i))){
#       results(i) = x(i);
#     } else {
#       results(i) = mu * mu / x(i);
#     }
#   }
#   return results;
# }")
#
#
# mu <- 0.07
# lambda <- 4
# x <- rinvgaussian_c(100000, mu, lambda)
# summary(x)
# quantile(x, probs = c(0.99))
# hist(x, nclass = 200, prob = TRUE)
# curve(exp(dinvgaussian(x, mu, lambda)), col = "red", add = TRUE)
