# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#

## define target distribution
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-2, 2), sd = 0.2, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}


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

### Test of the single kernel
niterations <- 10000
current_value <- rnorm(1, sd=10)
current_pdf <- target(current_value)
chain <- matrix(ncol=1, nrow=niterations)
chain[1,] <- current_value
for (t in 2:niterations){
  current_value <- single_kernel(current_value)
  chain[t,] <- current_value
}
qplot(x=1:niterations, y=chain, geom = "line") + ylab("X") + xlab("iteration")
hist(chain[1000:niterations], nclass = 100, prob = TRUE)
mean(chain > 1)

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

### test of the coupled kernel
niterations <- 10000
current_value1 <- rnorm(1, sd=10)
current_value2 <- rnorm(1, sd=10)
current_pdf1 <- target(current_value1)
current_pdf2 <- target(current_value2)
chain1 <- matrix(ncol=1, nrow=niterations)
chain2 <- matrix(ncol=1, nrow=niterations)
chain1[1,] <- current_value1
chain2[1,] <- current_value2
for (t in 2:niterations){
  current_value <- coupled_kernel(current_value1, current_value2)
  current_value1 <- current_value$chain_state1
  current_value2 <- current_value$chain_state2
  chain1[t,] <- current_value1
  chain2[t,] <- current_value2
}
it <- floor(seq(from = 1, to  = niterations, length.out = 1000))
qplot(x=it, y=chain1[it], geom = "line") + ylab("X") + xlab("iteration") + geom_line(aes(y = chain2[it]), colour = "red")

hist(chain1[1000:niterations], nclass = 100, prob = TRUE)
hist(chain2[1000:niterations], nclass = 100, prob = TRUE, add = TRUE, col = rgb(1,0,0,0.5))
curve(exp(target(x)), add = T)

## initial distribution
rinit <- function() rnorm(1)
res <- coupled_chains(single_kernel, coupled_kernel, rinit, K = 3)
res$meetingtime
res$iteration
H_bar(res, k=3, K = 3)
# #
qplot(x = 0:res$iteration, y = res$samples1[,1], geom = "line") +
  geom_line(aes(x = 0:(res$iteration-1), y = res$samples2[,1]), col = rgb(1,0,0,0.5)) +
  geom_line(aes(x = 1:res$iteration, y = res$samples2[,1]), col = rgb(1,0,0,1), linetype = 2)

### distribution of meeting times
nsamples <- 5000
meetingtime <- foreach(irep = 1:nsamples, .combine = c) %dorng% {
  c_chains <- coupled_chains(single_kernel, coupled_kernel, rinit, K = 1)
  c_chains$meetingtime
}
summary(meetingtime)
hist(meetingtime)
##
K <- 30
k <- 20
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
}
# histogram
component <- 1
find_breaks <- function(c_chains, component, nclass){
  all_samples <- unlist(sapply(c_chains, function(x) x$samples1[,component]))
  all_samples <- c(all_samples, unlist(sapply(c_chains, function(x) x$samples2[,component])))
  br <- hist(all_samples, plot=F, nclass = nclass)$breaks
  return(br)
}

nclass <- 20
breaks <- find_breaks(c_chains_, 1, nclass)
mids <- c()
for (i in 2:length(breaks)){
  mids <- c(mids, breaks[i-1] + (breaks[i] - breaks[i-1])/2)
}
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
