# load packages
library(unbiasedmcmc)
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)

target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

get_pb <- function(sd_proposal, initmean, initsd){
  # Markov kernel of the chain
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

  # Markov kernel of the coupled chain
  coupled_kernel <- function(state1, state2){
    chain_state1 <- state1$chain_state; current_pdf1 <- state1$current_pdf
    chain_state2 <- state2$chain_state; current_pdf2 <- state2$current_pdf
    proposal_value <- rnorm_max_coupling(chain_state1, chain_state2, sd_proposal, sd_proposal)
    proposal1 <- proposal_value$xy[1]
    proposal2 <- proposal_value$xy[2]
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- proposal_pdf1
    if (!proposal_value$identical){
      proposal_pdf2 <- target(proposal2)
    }
    logu <- log(runif(1))
    accept1 <- FALSE; accept2 <- FALSE
    if (is.finite(proposal_pdf1)) accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    if (is.finite(proposal_pdf2)) accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    if (accept1){
      chain_state1 <- proposal1
      current_pdf1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    identical_ <- ((proposal_value$identical) && (accept1) && (accept2))
    return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
                state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
                identical = identical_))
  }
  rinit <- function(){
    chain_state <- rnorm(1, initmean, initsd)
    current_pdf <- target(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  return(list(rinit = rinit, single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}

pb <- get_pb(3, initmean = 10, initsd = 10)

## test
# state1 <- pb$rinit()
# state2 <- pb$rinit()
# state1 <- pb$single_kernel(state1)
# pb$coupled_kernel(state1, state2)

# easy setting: good proposal
nsamples <- 1000
meetingtimes.easy <-  foreach(irep = 1:nsamples, .combine = c) %dorng% {
  sample_meetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit)$meetingtime
}
# hist(meetingtimes.easy)

m <- 10000
c_chains.easy <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
save(meetingtimes.easy, c_chains.easy, file = "bimodal.c_chains.easy.RData")


# intermediate setting: bad proposal, good init
pb <- get_pb(1, initmean = 10, initsd = 10)
meetingtimes.intermediate <-  foreach(irep = 1:nsamples, .combine = c) %dorng% {
  sample_meetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit)$meetingtime
}
# hist(meetingtimes.intermediate)

m <- 30000
c_chains.intermediate <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
save(meetingtimes.intermediate, c_chains.intermediate, file = "bimodal.c_chains.intermediate.RData")


# hard setting: bad proposal, bad init
pb <- get_pb(1, initmean = 10, initsd = 1)
meetingtimes.hard <-  foreach(irep = 1:nsamples, .combine = c) %dorng% {
  sample_meetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit)$meetingtime
}
# hist(meetingtimes.hard)

## coupled chains, first 1000 repeats
m <- 1000
nsamples <- 1000
c_chains.hard <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
save(meetingtimes.hard, c_chains.hard, file = "bimodal.c_chains.hard.RData")

## coupled chains, next 9000 repeats
m <- 1000
nsamples <- 9000
c_chains.hard.extra <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}

save(meetingtimes.hard, c_chains.hard, c_chains.hard.extra, file = "bimodal.c_chains.hard.RData")
