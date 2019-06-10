### This script plays with parallel tempering
### for a basic Ising model, with different values of the temperatures
### i.e. Gibbs sampler at each temperature and swap moves between chains
library(debiasedmcmc)
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
registerDoParallel(cores = detectCores()-2)

# size of the grid
size <- 32
# possible values of sum of neighbors
ss_ <- c(-4,-2,0,2,4)

nchains <- 16
# number of iterations
niterations <- 5e4
# history
history_sumstates <- matrix(0, ncol = nchains, nrow = niterations)
# probability of doing a swap move
proba_swapmove <- 0.01
# inverse temperatures
betas <- seq(from = 0.3, to = 0.55, length.out = nchains)
# precomputed probability for single-site flips
probas_ <-  sapply(betas, function(beta) exp(ss_*beta) / (exp(ss_*beta) + exp(-ss_*beta)))
# initialization
current_states <- ising_pt_rinit(nchains)
sumstates <- unlist(lapply(current_states, debiasedmcmc:::ising_sum_))
history_sumstates[1,] <- sumstates
nswap_attempts <- 0
nswap_accepts <- rep(0, nchains-1)
# MCMC loop
for (iteration in 2:niterations){
  if (iteration %% 10000 == 1) cat("iteration", iteration, "/", niterations, "\n")
  res_ <- ising_pt_single_kernel(current_states, sumstates, betas, probas_, proba_swapmove)
  current_states <- res_$chain_states
  sumstates <- res_$sumstates
  nswap_accepts <- nswap_accepts + res_$nswap_accepts
  nswap_attempts <- nswap_attempts + res_$nswap_attempts
  history_sumstates[iteration,] <- sumstates
}
save(nswap_accepts, nswap_attempts, history_sumstates, file = "ising.mcmc.RData")
load("ising.mcmc.RData")

niterations <- nrow(history_sumstates)
# swap acceptance in %
cat(paste0(round(100*nswap_accepts/nswap_attempts, 2), " %"), "\n")
par(mfrow = c(2,1))
# traceplots
matplot(apply(history_sumstates[2000:niterations,1:5], 2, function(x) cumsum(x)/(1:(niterations-2000+1))), type = "l", ylab = "sum of states at different temperatures")
# plot estimated average natural parameters
plot(betas, colMeans(history_sumstates[2000:niterations,]), type = "l")

par(mfrow = c(1,1))
matplot(history_sumstates[2000:1e4,1:5], type = "l")

library(coda)
vars <- sapply(1:nchains, function(i) spectrum0(history_sumstates[1e4:niterations,i])$spec)
sd_errors <- sqrt(vars/(niterations-1e4))
sd_errors
estimates <- colMeans(history_sumstates[1e4:niterations,])
ggplot(data.frame(betas = betas, estimates = estimates, sd_errors = sd_errors),
       aes(x = betas, y = estimates, ymin = estimates-2*sd_errors, ymax=estimates+2*sd_errors)) + geom_errorbar() +
  geom_line()
