
# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
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

coupled_kernel <- function(state1, state2){
  lambda1 <- state1$chain_state[1:ndata]
  beta1 <- state1$chain_state[ndata+1]
  lambda2 <- state2$chain_state[1:ndata]
  beta2 <- state2$chain_state[ndata+1]
  identical_components <- logical(length = ndata + 1)
  for (k in 1:ndata){
    x <- rgamma_coupled(alpha1 = alpha + s[k], alpha2 = alpha + s[k],
                        beta1 = beta1 + t[k], beta2 = beta2 + t[k])
    lambda1[k] <- x$xy[1]
    lambda2[k] <- x$xy[2]
    identical_components[k] <- x$identical
  }
  x <- rgamma_coupled(alpha1 = gamma + 10*alpha, alpha2 = gamma + 10*alpha,
                      beta1 = delta + sum(lambda1), beta2 = delta + sum(lambda2))
  beta1 <- x$xy[1]
  beta2 <- x$xy[2]
  identical_components[ndata+1] <- x$identical
  chain_state1 = c(lambda1, beta1); chain_state2 = c(lambda2, beta2)
  identical_ <- all(identical_components)
  return(list(state1 = list(chain_state = chain_state1), state2 = list(chain_state = chain_state2), identical = identical_))
}

k <- 7
m <- 10*7

nsamples <- 1000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(single_kernel, coupled_kernel, rinit, m=m)
}
save(nsamples, k, m, c_chains_, file = "pump.c_chain.RData")

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
table(meetingtime)
