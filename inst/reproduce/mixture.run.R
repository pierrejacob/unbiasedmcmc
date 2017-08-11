# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#

## target distribution
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

curve(exp(sapply(x, target)), from = -15, to = 15)

sd_proposal <- 3
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

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2){
  proposal_value <- rnorm_max_coupling(chain_state1, chain_state2, sd_proposal, sd_proposal)
  proposal1 <- proposal_value[1]
  proposal2 <- proposal_value[2]
  proposal_pdf1 <- target(proposal1)
  proposal_pdf2 <- target(proposal2)
  current_pdf1 <- target(chain_state1)
  current_pdf2 <- target(chain_state2)
  logu <- log(runif(1))
  accept1 <- FALSE
  accept2 <- FALSE
  if (is.finite(proposal_pdf1)){
    accept1 <- (logu < (proposal_pdf1 - current_pdf1))
  }
  if (is.finite(proposal_pdf2)){
    accept2 <- (logu < (proposal_pdf2 - current_pdf2))
  }
  if (accept1){
    chain_state1 <- proposal1
  }
  if (accept2){
    chain_state2 <- proposal2
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
}

rinit <- function() rnorm(1, mean = 10)
nsamples <- 10000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit)
}
save(c_chains_, file = "mixture.c_chains.RData")
load(file = "mixture.c_chains.RData")

##
K <- 200
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}
save(c_chains_, c_chains_continued_, file = "mixture.c_chains.RData")
load(file = "mixture.c_chains.RData")


sd_proposal <- 1
c_chains_2 <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit)
}
save(c_chains_, c_chains_continued_, c_chains_2, file = "mixture.c_chains.RData")
load(file = "mixture.c_chains.RData")
#
c_chains_continued_2 <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_2[[irep]], single_kernel, K = K)
}
save(c_chains_, c_chains_continued_, c_chains_2, c_chains_continued_2, file = "mixture.c_chains.RData")
load(file = "mixture.c_chains.RData")

