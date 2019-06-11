library(debiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

## this script verifies the validity of the H_bar function to compute
## estimators based on the run of coupled chains

target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-2, 2), sd = 0.2, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

curve(sapply(X = x, FUN = function(xp) exp(target(xp))), from = -5, to = 5)

# Markov kernel of the chain: MH with a Normal random walk proposal
Sigma_proposal <- diag(2^2, 1, 1)

# get MH kernels
kernels <- get_mh_kernels(target, Sigma_proposal)
single_kernel <- kernels$single_kernel
coupled_kernel <- kernels$coupled_kernel

rinit <- function(){
  x <- rnorm(1)
  return(list(chain_state = x, current_pdf = target(x)))
}

# kernels$single_kernel(rinit())
# kernels$coupled_kernel(rinit(), rinit())

### distribution of meeting times
nsamples <- 500
##
meetings_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit)
}

hist(sapply(meetings_, function(x) x$meetingtime))
k <- as.numeric(quantile(x = sapply(meetings_, function(x) x$meetingtime), probs = 0.95))
m <- 5 * k

c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m)
}

index_ <- which(sapply(c_chains_, function(x) x$meetingtime) > k)[1]

H_bar_test <- function(c_chains, h = function(x) x, k = 0, m = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if (m > maxiter){
    print("error: m has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  deltas <- matrix(0, nrow = maxiter+1, ncol = p) # p is number of columns is dim(image(h))
  deltas_term <- rep(0, p)
  for (t in 1:maxiter){ # t is as in the report, where the chains start at t=0
    deltas[t,] <- h(c_chains$samples1[t + 1,]) - h(c_chains$samples2[t,])
    if ((t >= (k+1)) && (t <= m)){
      deltas_term <- deltas_term + deltas[t,] * (t - k)
    }
  }
  deltas_term <- deltas_term + (m - k + 1) * (apply(X = deltas[(m+1):(maxiter+1),,drop=FALSE], 2, sum))
  h_of_chain <- apply(X = c_chains$samples1[(k+1):(m+1),,drop=F], MARGIN = 1, FUN = h)
  if (is.null(dim(h_of_chain))){
    h_of_chain <- matrix(h_of_chain, ncol = 1)
  } else {
    h_of_chain <- t(h_of_chain)
  }
  H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
  H_bar <- H_bar + deltas_term
  return(H_bar / (m - k + 1))
}

H_bar_test(c_chains_[[index_]], h = function(x) x, k = k, m = m)
H_bar(c_chains_[[index_]], h = function(x) x, k = k, m = m)

# Function that computes H_k for a given k
H_k <- function(c_chains, h = function(x) x, k = 0){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  p <- ncol(c_chains$samples1)
  deltas <- matrix(0, nrow = maxiter, ncol = p) # p is number of columns is dim(image(h))
  for (t in 1:maxiter){ # t is as in the report, where the chains start at t=0
    deltas[t,] <- h(c_chains$samples1[t + 1,]) - h(c_chains$samples2[t,])
  }
  H_k <- h(c_chains$samples1[k + 1,])
  if (k+1<=maxiter){
    H_k <- H_k + apply(X = deltas[((k+1):maxiter),,drop=FALSE], MARGIN = 2, FUN = sum)
  }
  return(H_k)
}

# function that computes H_bar by computing H_k for a range of values of k
H_bar_bruteforce <- function(c_chains, h = function(x) x, k = 0, m = 1){
  if (k >= m){
    print("error: k has to be < m")
    return(NULL)
  }
  H_bar <- H_k(c_chains, k = k)
  p <- length(H_bar)
  for (i in (k+1):m){
    H_bar <- H_bar + H_k(c_chains, k = i)
  }
  return(H_bar / (m-k+1))
}


H_bar_alternative <- function(c_chains, h = function(x) x, k = 0, m = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if (m > maxiter){
    print("error: m has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  h_of_chain <- apply(X = c_chains$samples1[(k+1):(m+1),,drop=F], MARGIN = 1, FUN = h)
  if (is.null(dim(h_of_chain))){
    h_of_chain <- matrix(h_of_chain, ncol = 1)
  } else {
    h_of_chain <- t(h_of_chain)
  }
  H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
  if (c_chains$meetingtime <= k + 1){
    # nothing else to add
  } else {
    deltas <- matrix(0, nrow = maxiter - k + 1, ncol = p)
    deltas_term <- rep(0, p)
    for (t in k:min(maxiter-1, c_chains$meetingtime-1)){ # t is as in the report, where the chains start at t=0
      coefficient <- min(t - k + 1, m - k + 1)
      delta_tp1 <- h(c_chains$samples1[t + 1 + 1,]) - h(c_chains$samples2[t+1,]) # the +1's are because R starts indexing at 1
      deltas_term <- deltas_term + coefficient * delta_tp1
    }
    H_bar <- H_bar + deltas_term
  }
  return(H_bar / (m - k + 1))
}


H_bar(c_chains_[[1]], h = function(x) x, k = k, m = m)
H_bar_alternative(c_chains_[[1]], h = function(x) x, k = k, m = m)
H_bar_bruteforce(c_chains_[[1]], k = k, m = m)

H_bar(c_chains_[[index_]], h = function(x) x, k = k, m = m)
H_bar_alternative(c_chains_[[index_]], h = function(x) x, k = k, m = m)
H_bar_bruteforce(c_chains_[[index_]], k = k, m = m)

test_H_bar <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
  H_bar(c_chains_[[irep]], h = function(x) x, k = k, m = m)
}

test_H_bar_alternative <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
  H_bar_alternative(c_chains_[[irep]], h = function(x) x, k = k, m = m)
}

test_H_bar_bruteforce <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
  H_bar_bruteforce(c_chains_[[irep]], h = function(x) x, k = k, m = m)
}

summary(as.numeric(abs(test_H_bar - test_H_bar_alternative)))
summary(as.numeric(abs(test_H_bar - test_H_bar_bruteforce)))

mean(test_H_bar)
var(test_H_bar)
