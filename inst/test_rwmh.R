# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#

target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-2, 2), sd = 0.2, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

# curve(target(x), from = -5, to = 5)
# Markov kernel of the chain
single_kernel <- function(chain_state){
  proposal_value <- rnorm(1, mean=chain_state, sd=2)
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
Sigma <- diag(2^2, 1, 1)
coupled_kernel <- function(chain_state1, chain_state2){
  proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma, Sigma)
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

rinit <- function() rnorm(1)
nsamples <- 10000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)
##
K <- 100
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}

summary(sapply(c_chains_continued_, function(x) x$meetingtime))
summary(sapply(c_chains_, function(x) x$iteration))

k <- 20
# histogram
histogram <- histogram_c_chains(c_chains_continued_, 1, k, K, nclass = 100)
g <- plot_histogram(histogram)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red")
g

