# load packages
library(debiasedmcmc)
rm(list = ls())
set.seed(21)

#
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
#
get_pb <- function(sd_proposal, initmean, initsd){
  # Markov kernel of the chain
  single_kernel <- function(chain_state){
    proposal_value <- rnorm(1, mean=chain_state, sd=sd_proposal)
    proposal_pdf <- target(proposal_value)
    current_pdf <- target(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(proposal_value)
    } else {
      return(chain_state)
    }
  }
  rinit <- function() rnorm(1, initmean, initsd)
  return(list(rinit = rinit, single_kernel = single_kernel))
}
#

# test function
testfunction <- function(x) (x > 3)
## Modify nmcmc
nmcmc <- 100000
burnin <- 10000
# easy setting: good proposal, but bad init
pb <- get_pb(3, initmean = 10, initsd = 10)
current <- pb$rinit()
chain <- rep(0, nmcmc)
for (imcmc in 1:nmcmc){
  current <- pb$single_kernel(current)
  chain[imcmc] <- current
}


library(coda)
mcmcvar.easy <- spectrum0(sapply(chain[burnin:nmcmc], testfunction))$spec
mcmcvar.easy
mean(sapply(chain[burnin:nmcmc], testfunction))
save(nmcmc, mcmcvar.easy, file = "mixture.mcmc.RData")


# intermediate setting: bad proposal, good init
pb <- get_pb(1, initmean = 10, initsd = 10)
current <- pb$rinit()
chain <- rep(0, nmcmc)
for (imcmc in 1:nmcmc){
  current <- pb$single_kernel(current)
  chain[imcmc] <- current
}
library(coda)
mcmcvar.intermediate <- spectrum0(sapply(chain[burnin:nmcmc], testfunction))$spec
save(nmcmc, mcmcvar.easy, mcmcvar.intermediate, file = "mixture.mcmc.RData")


