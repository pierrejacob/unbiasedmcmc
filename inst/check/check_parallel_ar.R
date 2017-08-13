# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
library(doRNG)


# logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
# kernel <- function(state) 0.5*state + rnorm(1, 0, sqrt(3/4))
# ckernel <- function(state1, state2){
#   states <- rnorm_max_coupling(.5*state1, .5*state2, sqrt(3/4), sqrt(3/4))
#   return(list(chain_state1 = states[1], chain_state2 = states[2]))
# }
#
# rinit <- function() rnorm(1, 10, 1)
library(parallel)
cl <- makeCluster(12, type = "FORK")
registerDoParallel(cl)
# clusterSetRNGStream(cl, 12)
nsamples <- 50000
c_chains_ <- foreach(icount(nsamples)) %dopar% {
  logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
  kernel <- function(state) 0.5*state + rnorm(1, 0, sqrt(3/4))
  ckernel <- function(state1, state2){
    states <- rnorm_max_coupling(.5*state1, .5*state2, sqrt(3/4), sqrt(3/4))
    return(list(chain_state1 = states[1], chain_state2 = states[2]))
  }
  rinit <- function() rnorm(1, 10, 1)
  coupled_chains(kernel, ckernel, rinit)
}
stopCluster(cl)

meetingtimes <- sapply(c_chains_, function(x) x$meetingtime)
print(summary(meetingtimes))
k <- as.numeric(quantile(meetingtimes, probs = 0.95))
K <- 10*k
# clusterExport(cl, "K")
c_chains_ <- foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(kernel, ckernel, rinit, K = K)
}

histogram1 <- histogram_c_chains(c_chains_, 1, k, k, nclass = 30)
plot_histogram(histogram1)


# stopCluster(cl)
