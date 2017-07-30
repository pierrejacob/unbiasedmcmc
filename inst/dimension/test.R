# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
setwd("~/Dropbox/PolyaGammaResults/dimension/")

## try coupling a MH for multivariate Normal directly
dimension <- 100
target <- function(x) fast_dmvnorm(matrix(x, nrow = 1), rep(0, dimension), diag(1, dimension, dimension))

# curve(target(x), from = -5, to = 5)
Sigma_proposal <- diag(1/dimension, dimension, dimension)
# Markov kernel of the chain
single_kernel <- function(chain_state){
  proposal_value <- fast_rmvnorm(1, chain_state, Sigma_proposal)
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
  proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
  proposal1 <- proposal_value[,1]
  proposal2 <- proposal_value[,2]
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

rinit <- function() fast_rmvnorm(1, rep(0, dimension), Sigma_proposal)
nsamples <- 10
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit, max_iterations = 1000)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
# hist(meetingtime)

# K <- 200
# c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
#   continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
# }
#
# k <- 50
# histogram <- histogram_c_chains(c_chains_continued_, 1, k, K, nclass = 100)
# g <- plot_histogram(histogram, with_bar = T)
# g <- g + stat_function(fun = function(x) dnorm(x), colour = "red", alpha = 1)
# g


iterate_per_component <- 1
## now try the MH within Gibbs way
gibbs_update <- function(component){
  single_step_ <- function(chain_state){
    increment <- rnorm(1)
    proposal_value <- chain_state
    proposal_value[component] <- proposal_value[component] + increment
    proposal_pdf <- target(proposal_value)
    current_pdf <- target(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(proposal_value)
    } else {
      return(chain_state)
    }
  }
  return(single_step_)
}

single_kernel <- function(chain_state){
  for (component in 1:dimension){
    step_ <- gibbs_update(component)
    for (iterate in 1:iterate_per_component){
      chain_state <- step_(chain_state)
    }
  }
  return(chain_state)
}
# niterations <- 1000
# chains <- matrix(nrow = niterations, ncol = dimension)
# state <- rinit()
# for (iteration in 1:niterations){
#   state <- single_kernel(state)
#   chains[iteration,] <- state
# }
# hist(chains[,1], prob = TRUE)
# curve(dnorm(x), add = T, col = "red")

coupled_gibbs_update <- function(component){
  coupled_kernel_ <- function(chain_state1, chain_state2){
    proposal_value <- gaussian_max_coupling(chain_state1[component], chain_state2[component], diag(1,1,1), diag(1,1,1))
    proposal1 <- chain_state1
    proposal1[component] <- proposal_value[,1]
    proposal2 <- chain_state2
    proposal2[component] <- proposal_value[,2]
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
  return(coupled_kernel_)
}

coupled_kernel <- function(chain_state1, chain_state2){
  for (component in 1:dimension){
    step_ <- coupled_gibbs_update(component)
    for (iterate in 1:iterate_per_component){
      chain_states <- step_(chain_state1, chain_state2)
      chain_state1 <- chain_states$chain_state1
      chain_state2 <- chain_states$chain_state2
    }
  }
  return(chain_states)
}

nsamples <- 100
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit, max_iterations = 1000)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)

K <- 100
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}

k <- 40
histogram <- histogram_c_chains(c_chains_continued_, 1, k, K, nclass = 100)
g <- plot_histogram(histogram, with_bar = T)
g <- g + stat_function(fun = function(x) dnorm(x), colour = "red", alpha = 1)
g


