
# load packages
library(unbiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())

#
##  This example is about failures of nuclear pumps. It's classic (e.g. Example 10.17 in Robert & Casella Monte Carlo Statistical Methods)
## It's used as an example in Murdoch and Green's perfect samplers paper and also Reutter and Johnson 1995
## about using coupled chains to monitor MCMC convergence

# The data:
# number of failures
s <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
# times
t <- c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
# the model says s_k ~ Poisson(lambda_k * t_k), for k = 1,...,10
ndata <- 10
# and lambda_k ~ Gamma(alpha,beta), beta ~ Gamma(gamma, delta)
alpha <- 1.802
gamma <- 0.01
delta <- 1
# full conditionasl:
# lambda_k given rest: Gamma(alpha + s_k, beta + t_k)
# beta given rest: Gamma(gamma + 10*alpha, delta + sum_{k=1}^10 lambda_k)

single_kernel <- function(state){
  lambda <- state$chain_state[1:ndata]
  beta <- state$chain_state[ndata+1]
  for (k in 1:ndata){
    lambda[k] <- rgamma(1, shape = alpha + s[k], rate = beta + t[k])
  }
  beta <- rgamma(1, shape = gamma + 10*alpha, rate = delta + sum(lambda))
  return(list(chain_state = c(lambda, beta)))
}

rinit <- function(){
  return(list(chain_state = rep(1, ndata+1)))
}


niterations <- 5e5
chain <- matrix(nrow = niterations, ncol = ndata+1)
current <- rinit()
chain[1,] <- current$chain_state
for (iteration in 2:niterations){
  current <- single_kernel(current)
  chain[iteration,] <- current$chain_state
}
save(niterations, chain, file = "pump.mcmc.RData")
# load("pump.mcmc.RData")
# hist(chain[,ndata+1])
