### This script plays with Gibbs sampler (i.e. single site updates)
### for a basic Ising model, with different values of the temperatures
library(debiasedmcmc)
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
registerDoParallel(cores = detectCores()-2)


## check rcouplbern2
# debiasedmcmc:::rcouplbern2(0.6, 0.7)
# xs <- foreach(i = 1:10000, .combine = rbind) %dorng% {
#   debiasedmcmc:::rcouplbern2(0.6, 0.8)
# }
# table(xs[,1])
# table(xs[,2])
# mean(apply(xs, 1, function(x) x[1] == x[2]))

# size of the grid
size <- 32
# possible values of sum of neighbors
ss_ <- c(-4,-2,0,2,4)
# inverse temperature
beta <- 0.35
# precomputed probability for single-site flips
proba_ <-  exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta))
# initialization
chain_state <- ising_rinit()
# MCMC loop
niterations <- 100
# store history
history_states <- array(data = 0, dim = c(niterations, size, size))
history_states[1,,] <- chain_state
#
for (iter in 2:niterations){
  chain_state <- ising_single_kernel(chain_state, proba_)
  history_states[iter,,] <- chain_state
}
# (peek)
# image(history_states[niterations,,])
#

# this computes the evolution of the sum of {x_i x_j} over neighbors,
# also called the "natural statistic" in Geyer 1991
sumstates <- apply(X = history_states, MARGIN = 1, FUN = debiasedmcmc:::ising_sum_)
plot(sumstates, type = "l")

# Now let's run Gibbs chains at different betas
run_gibbs <- function(niterations, beta){
  proba_ <-  exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta))
  chain_state <- ising_rinit()
  sumstates <- rep(0, niterations)
  sumstates[1] <- debiasedmcmc:::ising_sum_(chain_state)
  for (iter in 2:niterations){
    chain_state <- ising_single_kernel(chain_state, proba_)
    sumstates[iter] <- debiasedmcmc:::ising_sum_(chain_state)
  }
  return(sumstates)
}

niterations <- 5000
betas <- seq(from = 0.35, to = 0.50, length.out = 20)
sumstates <- matrix(ncol = length(betas), nrow = niterations)
for (ibeta in seq_along(betas)){
  sumstates[,ibeta] <- run_gibbs(niterations, betas[ibeta])
}

#
matplot(sumstates[,1], type = "l")
matplot(sumstates[,10], type = "l")
matplot(sumstates[,20], type = "l")

plot(betas, colMeans(sumstates[2000:niterations,]), type = "l")

