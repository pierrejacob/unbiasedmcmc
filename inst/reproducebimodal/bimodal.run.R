# load packages
library(debiasedmcmc)
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())


#
samplemeetingtime <- function(single_kernel, coupled_kernel, rinit, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  chain_state1 <- single_kernel(chain_state1)
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      chain_state1 <- single_kernel(chain_state1)
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    # stop after max(m, tau) steps
    if (iter >= meetingtime){
      finished <- TRUE
    }
  }
  return(list(meetingtime = meetingtime, iteration = iter, finished = finished))
}

## hard target distribution

target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

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
  rinit <- function() rnorm(1, initmean, initsd)
  return(list(rinit = rinit, single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}


# easy setting: good proposal
pb <- get_pb(3, initmean = 10, initsd = 10)
nsamples <- 1000
meetingtimes.easy <-  foreach(irep = 1:nsamples, .combine = c) %dorng% {
  samplemeetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit)$meetingtime
}
hist(meetingtimes.easy)

m <- 10000
nsamples <- 1000
c_chains.easy <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
save(meetingtimes.easy, c_chains.easy, file = "mixture.c_chains.easy.RData")


# intermediate setting: bad proposal, good init
pb <- get_pb(1, initmean = 10, initsd = 10)
nsamples <- 1000
meetingtimes.intermediate <-  foreach(irep = 1:nsamples, .combine = c) %dorng% {
  samplemeetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit)$meetingtime
}
hist(meetingtimes.intermediate)

m <- 30000
nsamples <- 1000
c_chains.intermediate <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
save(meetingtimes.intermediate, c_chains.intermediate, file = "mixture.c_chains.intermediate.RData")


# hard setting: bad proposal, bad init
pb <- get_pb(1, initmean = 10, initsd = 1)
nsamples <- 1000
meetingtimes.hard <-  foreach(irep = 1:nsamples, .combine = c) %dorng% {
  samplemeetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit)$meetingtime
}
hist(meetingtimes.hard)

m <- 1000
nsamples <- 1000
c_chains.hard <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
sapply(c_chains.hard, function(x) x$meetingtime) %>% hist
save(meetingtimes.hard, c_chains.hard, file = "mixture.c_chains.hard.RData")

m <- 1000
nsamples <- 9000
c_chains.hard.extra <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit, m = m)
}
sapply(c_chains.hard.extra, function(x) x$meetingtime) %>% hist
save(meetingtimes.hard, c_chains.hard, c_chains.hard.extra, file = "mixture.c_chains.hard.RData")



### MCMC
