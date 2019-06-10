# load packages
library(debiasedmcmc)
library(coda)
#
rm(list = ls())
set.seed(1)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores() - 1)

data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)


mcmc_blasso <- function(nmcmc, burnin, lambda){
  pb <- get_blasso(Y, X, lambda)
  state <- pb$rinit()
  chain <- matrix(nrow = nmcmc, ncol = length(state$chain_state))
  for (imcmc in 1:nmcmc){
    state <- pb$gibbs_kernel(state)
    chain[imcmc,] <- state$chain_state
  }
  return(chain)
}


## Modify nmcmc
nmcmc <- 5000
burnin <- floor(nmcmc / 10)

# result <- mcmc_blasso(nmcmc, burnin, .1)
# # matplot(result[,1:10], type = "l")
# postmeans <- colMeans(result[burnin:nmcmc,])
# ess <- effectiveSize(result[burnin:nmcmc,1:p])

lambdas <- 10^(seq(from = -2, to = 3, length.out = 25))
df <- foreach (ilambda = 1:length(lambdas), .combine = rbind) %dorng% {
  lambda <- lambdas[ilambda]
  print(lambda)
  result <- mcmc_blasso(nmcmc, burnin, lambda)
  postmeans <- colMeans(result[burnin:nmcmc,1:p])
  ess <- effectiveSize(result[burnin:nmcmc,1:p])

  data.frame(ilambda = rep(ilambda, p), lambda = rep(lambda, p), component = 1:p,
                   ess = ess, postmeans = postmeans)
}
save(df, nmcmc, burnin, lambdas, file = "bayesianlasso.mcmc.RData")


