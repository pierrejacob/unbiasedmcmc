library(unbiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
# remove all
rm(list = ls())
# set RNG seed
set.seed(11)

logtarget <- function(x) dnorm(x, mean = 0, sd = 1, log = TRUE)
rinit <- function(){
  chain_state <- rnorm(1, 1, 1)
  current_pdf <- logtarget(chain_state)
  return(list(chain_state = chain_state, current_pdf = current_pdf))
}
sd_proposal <- 1
MH_kernel <- function(state){
  chain_state <- state$chain_state
  current_pdf <- state$current_pdf
  proposal <- rnorm(1, chain_state, sd_proposal)
  proposal_pdf <- logtarget(proposal)
  if (log(runif(1)) < (proposal_pdf - current_pdf)){
    return(list(chain_state = proposal, current_pdf = proposal_pdf))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
}

niterations <- 10000
chain <- rep(0, niterations)
state <- rinit()
for (i in 1:niterations){
  state <- MH_kernel(state)
  chain[i] <- state$chain_state
}
hist(chain, prob = TRUE, nclass = 40, main = "")
curve(exp(logtarget(x)), add = TRUE, col = "red")

coupledMH_kernel <- function(state1, state2){
  chain_state1 <- state1$chain_state;  current_pdf1 <- state1$current_pdf
  chain_state2 <- state2$chain_state;  current_pdf2 <- state2$current_pdf
  # proposal from a maximal coupling
  proposal <- rnorm_max_coupling(chain_state1, chain_state2, sd_proposal, sd_proposal)
  proposal_pdf1 <- logtarget(proposal$xy[1])
  # only compute target pdf on 2nd proposal if it is not identical to 1st proposal
  proposal_pdf2 <- proposal_pdf1
  if (!proposal$identical){
    proposal_pdf2 <- logtarget(proposal$xy[2])
  }
  logu <- log(runif(1))
  accept1 <- FALSE; accept2 <- FALSE
  if (is.finite(proposal_pdf1)){
    if (logu < (proposal_pdf1 - current_pdf1)){
      accept1 <- TRUE
      chain_state1 <- proposal$xy[1]; current_pdf1 <- proposal_pdf1
    }
  }
  if (is.finite(proposal_pdf2)){
    if(logu < (proposal_pdf2 - current_pdf2)){
      accept2 <- TRUE
      chain_state2 <- proposal$xy[2]; current_pdf2 <- proposal_pdf2
    }
  }
  identical_ <- (proposal$identical) && (accept1) && accept2
  return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
              state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
              identical = identical_))
}

nsamples <- 10000
lag <- 10
k <- 10
m <- 50
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(MH_kernel, coupledMH_kernel, rinit, lag = lag, m = m)
}
names(c_chains_[[1]])
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
hist(meetingtime, breaks = 1:max(meetingtime), prob = TRUE, main = "", xlab = "meeting times")
mean_cost <- mean(sapply(c_chains_, function(x) x$cost))
uestimators <-  sapply(c_chains_, function(x) H_bar(x, h = function(x) x^2, k = k, m = m))

uestimators2 <-  foreach(irep = 1:nsamples) %dorng% {
  sample_unbiasedestimator(MH_kernel, coupledMH_kernel, rinit, h = function(x) x^2, k = k, m = m, lag = lag)
}

mean(uestimators)
mean(sapply(uestimators2, function(x) x$uestimator))
var(uestimators)
var(sapply(uestimators2, function(x) x$uestimator))

hist(uestimators)
hist(sapply(uestimators2, function(x) x$uestimator))

mean_cost
mean(sapply(uestimators2, function(x) x$cost))

# ## cost: lag calls to kernel, tau - lag calls to coupled kernel, then max(0, m - tau) calls to kernel again (if meeting < m)
# cost_of_coupled_chain <- function(c_chains){
#   lag <- dim(c_chains$samples1)[1] - dim(c_chains$samples2)[1]
#   cost <- lag + 2*(c_chains$meetingtime - lag) + max(0, c_chains$iteration - c_chains$meetingtime)
#   return(cost)
# }
#
# costs_ <- sapply(c_chains_, function(x) cost_of_coupled_chain(x))

inefficiency <- mean_cost * var(uestimators)
cat("inefficiency:", inefficiency, ", cost:", mean_cost)


