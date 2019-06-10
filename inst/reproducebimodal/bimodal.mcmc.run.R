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
  single_kernel <- function(state){
    chain_state <- state$chain_state
    current_pdf <- state$current_pdf
    proposal_value <- rnorm(1, mean=chain_state, sd=sd_proposal)
    proposal_pdf <- target(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
  }
  rinit <- function(){
    chain_state <- rnorm(1, initmean, initsd)
    current_pdf <- target(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }

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
  chain[imcmc] <- current$chain_state
}


library(coda)
mcmcvar.easy <- spectrum0(sapply(chain[burnin:nmcmc], testfunction))$spec
mcmcvar.easy
mean(sapply(chain[burnin:nmcmc], testfunction))
save(nmcmc, mcmcvar.easy, file = "bimodal.mcmc.RData")


# intermediate setting: bad proposal, good init
pb <- get_pb(1, initmean = 10, initsd = 10)
current <- pb$rinit()
chain <- rep(0, nmcmc)
for (imcmc in 1:nmcmc){
  current <- pb$single_kernel(current)
  chain[imcmc] <- current$chain_state
}
library(coda)
mcmcvar.intermediate <- spectrum0(sapply(chain[burnin:nmcmc], testfunction))$spec
save(nmcmc, mcmcvar.easy, mcmcvar.intermediate, file = "bimodal.mcmc.RData")


