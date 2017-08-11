# load packages
library(debiasedmcmc)
library(ggthemes)
library(lars)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
# registerDoParallel(cores = 1)
setwd("~/Dropbox/PolyaGammaResults/reproduce/")

data(diabetes)
X <- scale(diabetes$x)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

load(file = "diabetes.cchains.RData")

lambdas
head(df)

nmcmc <- 2000
burnin <- 1000
mcmc.df <- foreach (ilambda = 1:length(lambdas), .combine = rbind) %dorng% {
  lambda <- lambdas[ilambda]
  pb <- get_blasso(Y, X, lambda)
  chain_start <- pb$rinit()
  chain <- matrix(nrow = nmcmc, ncol = length(chain_start))
  chain[1,] <- chain_start
  for (iteration in 2:nmcmc){
    chain[iteration,] <- pb$gibbs_kernel(chain[iteration-1,])
  }
  mcmc_mean <- colMeans(chain[burnin:nmcmc,1:p])
  data.frame(ilambda = rep(ilambda, p), lambda = rep(lambda, p), component = 1:p,
             mcmc_mean = matrix(mcmc_mean, nrow = p))
}

