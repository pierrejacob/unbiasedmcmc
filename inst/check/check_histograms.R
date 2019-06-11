library(debiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

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
nsamples <- 1000
meetingtime <- foreach(irep = 1:nsamples, .combine = c) %dorng% {
  c_chains <- coupled_chains(single_kernel, coupled_kernel, rinit, K = 1)
  c_chains$meetingtime
}
summary(meetingtime)
hist(meetingtime)
##
K <- 50
k <- 20
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
}

index_ <- which(sapply(c_chains_, function(x) x$meetingtime) > 30)[1]
c_chains_[[index_]]
H_bar(c_chains_[[index_]], h = function(x) x, k = k, K = K)
#
summary(sapply(c_chains_, function(x) x$meetingtime))
# histogram
component <- 1

find_breaks <- function(c_chains, component, nclass){
  all_samples <- unlist(sapply(c_chains, function(x) x$samples1[,component]))
  all_samples <- c(all_samples, unlist(sapply(c_chains, function(x) x$samples2[,component])))
  br <- hist(all_samples, plot=F, nclass = nclass)$breaks
  return(br)
}

nclass <- 30
breaks <- find_breaks(c_chains_, 1, nclass)

mids <- c()
for (i in 2:length(breaks)){
  mids <- c(mids, breaks[i-1] + (breaks[i] - breaks[i-1])/2)
}
mids
width <- diff(breaks)[1]
# ## compute histogram
prop <- rep(0, length(breaks)-1)
sd_prop <- rep(0, length(breaks)-1)
for (ibreak in 2:length(breaks)){
  print(ibreak)
  estimators <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    H_bar(c_chains_[[irep]], h = function(x) (x[1]> breaks[ibreak-1] && x[1]<breaks[ibreak]), k = k, K = K)
  }
  prop[ibreak-1] <- mean(estimators)
  sd_prop[ibreak-1] <- sd(estimators) / sqrt(nsamples)
}

plot(x = mids, y = (prop+2*sd_prop) / width, type = "l")
lines(x = mids, y = prop / width, lty = 2)
lines(x = mids, y = (prop-2*sd_prop) / width)
curve(sapply(x, function(x_)exp(target(x_))), add = T, col = "red")

integrate(function(x) sapply(x, function(x_)exp(target(x_))), lower = -7, upper = 7)

ibreak <- 27
irep <- which(sapply(c_chains_, function(x) x$meetingtime > k))[1]
H_bar(c_chains_[[irep]], h = function(x) (x[1]> breaks[ibreak-1] && x[1]<breaks[ibreak]), k = k, K = K)
c_chains_[[irep]]$samples1
c_chains_[[1]]$samples2

cppFunction('double histogram_bin(List c_chains, int component, double lower, double upper, int k, int K){
  int meetingtime = c_chains["meetingtime"];
  int iteration = c_chains["iteration"];
  NumericMatrix samples1 = c_chains["samples1"];
  NumericMatrix samples2 = c_chains["samples2"];
  double estimator = 0;
  for (int isample = k; isample <= K; isample ++){
    //std::cerr << isample << ":" << samples1(isample,component-1) << std::endl;
    if (samples1(isample,component-1) > lower && samples1(isample,component-1) < upper){
      estimator += 1;
    }
  }
  if (meetingtime > k + 1){
    double coefficient = 0.;
    double increment = 0;
    for (int isample = k; isample <= std::min(iteration-1, meetingtime-1); isample ++){
      increment = 0;
      coefficient = std::min(isample - k + 1, K - k + 1);
      if (samples1(isample+1,component-1) > lower && samples1(isample+1,component-1) < upper){
        increment += coefficient;
      }
      if (samples2(isample,component-1) > lower && samples2(isample,component-1) < upper){
        increment -= coefficient;
      }
      estimator += increment;
    }
  }
  return estimator / (K - k + 1);
}')

histogram_bin(c_chains_[[irep]], 1, breaks[ibreak-1], breaks[ibreak], k, K)
H_bar(c_chains_[[irep]], h = function(x) (x[1]> breaks[ibreak-1] && x[1]<breaks[ibreak]), k = k, K = K)


prop2 <- rep(0, length(breaks)-1)
sd_prop2 <- rep(0, length(breaks)-1)
for (ibreak in 2:length(breaks)){
  print(ibreak)
  estimators <- foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    histogram_bin(c_chains_[[irep]], 1, breaks[ibreak-1], breaks[ibreak], k, K)
  }
  prop2[ibreak-1] <- mean(estimators)
  sd_prop2[ibreak-1] <- sd(estimators) / sqrt(nsamples)
}

prop - prop2
sd_prop - sd_prop2
