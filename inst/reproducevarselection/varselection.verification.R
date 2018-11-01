library(debiasedmcmc)
rm(list = ls())
set.seed(111)

# correlated design
n <- 500
p <- 1000
SNR <- 2
s_star <- 10
s0 <- 100
sigma0 <- 1
beta_star <- SNR * sqrt(sigma0^2 * log(p) / n) * c(2,-3,2,2,-3,3,-2,3,-2,3, rep(0, p-10))
#
# covariance <- matrix(0, nrow = p, ncol = p)
# for (i in 1:p){
#   for (j in 1:p){
#     covariance[i,j] <- exp(-abs(i-j))
#   }
# }
#
# X <- fast_rmvnorm(n, mean = rep(0, p), covariance)
# # X <- matrix(rnorm(n * p), nrow = n, ncol = p) # fast_rmvnorm_chol(n, rep(0, p), diag(1, p, p))
# X <- scale(X)
# Y <- X %*% matrix(beta_star, ncol = 1) + rnorm(n, 0, sigma0)
# Y <- scale(Y)
#
# save(X, Y, beta_star, file = paste0("varselection.data.verification.RData"))
load(paste0("varselection.data.verification.RData"))

## Now produce meeting times

g <- p^3

s0 <- 100
kappa <- 2
proportion_singleflip <- 0.5

vs <- get_variableselection(Y,X,g,kappa,s0,proportion_singleflip)
prior <- vs$prior
marginal_likelihood <- vs$marginal_likelihood
rinit <- vs$rinit
single_kernel <- vs$single_kernel
coupled_kernel <- vs$coupled_kernel
unbiasedestimator <- vs$unbiasedestimator
# Unbiased MCMC
nrep <- 100

library(doRNG)
library(doParallel)
registerDoParallel(cores = detectCores()-2)
ues <- foreach(i = 1:nrep) %dorng% {
  unbiasedestimator(single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 0)
}
meetings <- sapply(ues, function(x) x$meetingtime)
summary(meetings)



