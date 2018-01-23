# load packages
library(debiasedmcmc)
library(ggthemes)
library(coda)

#
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores() - 1)

data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)


mcmc_blasso <- function(nmcmc, burnin, lambda){
  pb <- get_blasso(Y, X, lambda)
  state <- pb$rinit()
  states <- matrix(nrow = nmcmc, ncol = length(state))
  for (imcmc in 1:nmcmc){
    state <- pb$gibbs_kernel(state)
    states[imcmc,] <- state
  }
  # ess <- try(effectiveSize(states[burnin:nmcmc,]))
  # if (inherits(ess, "try-error")){
  #   return(list(ess = rep(NA, 2*p+1), postmeans = rep(NA, 2*p+1)))
  # } else {
  #   postmeans <- colMeans(states[burnin:nmcmc,])
  #   return(list(ess = ess, postmeans = postmeans))
  # }
  return(states)
}

# lambda <- 1000
nmcmc <- 50000
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
save(df, lambdas, file = "bayesianlasso.mcmc.RData")
load("bayesianlasso.mcmc.RData")

# save(df, lambdas, nmcmc, burnin, file = "diabetes.p64.gibbsperformance.RData")
# df <- data.frame()
# for(ilambda in 1:length(lambdas)){
#   lambda <- lambdas[ilambda]
#   print(lambda)
#   result <- mcmc_blasso(nmcmc, burnin, lambda)
#   l <- data.frame(ilambda = rep(ilambda, 2*p+1), lambda = rep(lambda, 2*p+1), component = 1:(2*p+1),
#                    ess = result$ess, postmeans = result$postmeans)
#   df <- rbind(df, l)
#   save(df, lambdas, nmcmc, burnin, file = "bayesianlasso.mcmc.RData")
# }
# save(df, lambdas, nmcmc, burnin, file = "diabetes.p64.gibbsperformance.RData")
# load(file = "diabetes.p64.gibbsperformance.RData")


head(df)

ggplot(df, aes(x = lambda, y = postmeans)) + geom_point() + scale_x_log10() + geom_line(aes(group = component))

ggplot(df, aes(x = lambda, y = ess / (nmcmc - burnin + 1))) + geom_point() + scale_x_log10()

