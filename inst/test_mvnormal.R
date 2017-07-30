
# This code sets up a version of our estimator for drawing from a specified MVN distribution
# using (coupling) random walk MH. When the chains are far apart we use an optimal transport
# coupling, but when they are close an optimal coupling is used. The latter part of this code
# tests the chain continuation and histogram functions.

# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
dimension <- 4
mean_target <- rnorm(dimension)
Sigma_target <- diag(1, dimension, dimension)
for (i in 1:(dimension-1)){
  Sigma_target[i,i+1] <- Sigma_target[i+1,i] <- -0.3
}
target <- function(x){
  return(fast_dmvnorm(matrix(x, nrow = 1), mean_target, Sigma_target))
}

# curve(target(x), from = -5, to = 5)
# Markov kernel of the chain
# sd_proposal <- .5
single_kernel <- function(chain_state){
  proposal_value <- chain_state + fast_rmvnorm(1, mean_target, Sigma_target)[1,]
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
Sigma <- Sigma_target
Sigma_chol <- chol(Sigma)
Sigma_chol_inv <- solve(chol(Sigma))
coupled_kernel <- function(chain_state1, chain_state2){
  distance_ <- mean((chain_state1 - chain_state2)^2)
  if (distance_ > 5){
    proposal_value <- gaussian_opt_transport(1, chain_state1, chain_state2, Sigma_chol, Sigma_chol, Sigma_chol_inv, Sigma_chol_inv)[[1]]
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
  } else {
    proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma, Sigma)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
  }
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

rinit <- function() fast_rmvnorm(1, mean_target, Sigma_target)[1,]

### test of the coupled kernel
niterations <- 10000
current_value1 <- rinit()
current_value2 <- rinit()
current_pdf1 <- target(current_value1)
current_pdf2 <- target(current_value2)
chain1 <- matrix(ncol=dimension, nrow=niterations)
chain2 <- matrix(ncol=dimension, nrow=niterations)
chain1[1,] <- current_value1
chain2[1,] <- current_value2
for (t in 2:niterations){
  current_value <- coupled_kernel(current_value1, current_value2)
  current_value1 <- current_value$chain_state1
  current_value2 <- current_value$chain_state2
  chain1[t,] <- current_value1
  chain2[t,] <- current_value2
}

distances_ <- sapply(1:niterations, function(t) mean((chain1[t,] - chain2[t,])^2))
mean(diff(chain1[,1]) > 1e-10)
plot(1:niterations, distances_, type = "l")

it <- floor(seq(from = 1, to  = niterations, length.out = 1000))
qplot(x=it, y=chain1[it,1], geom = "line") + ylab("X") + xlab("iteration") + geom_line(aes(y = chain2[it,1]), colour = "red")
qplot(x=it, y=chain1[it,2], geom = "line") + ylab("X") + xlab("iteration") + geom_line(aes(y = chain2[it,2]), colour = "red")
acf(chain1[,1])
tail(chain1)
tail(chain2)

#
#
nsamples <- 1000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)

# ##
K <- 1000
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}
#
# summary(sapply(c_chains_continued_, function(x) x$meetingtime))
# summary(sapply(c_chains_, function(x) x$iteration))
#

k <- 700
# histogram
histogram <- histogram_c_chains(c_chains_continued_, 1, k, K, nclass = 100)
g <- plot_histogram(histogram)
g <- g + stat_function(fun = function(x) dnorm(x, mean = mean_target[1], sd = sqrt(Sigma_target[1,1])), colour = "red")
g

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], k = 700, K = 1000)
}

for (component in 1:dimension){
  estimators <- sapply(mean_estimators, function(x) x[component])
  cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/length(estimators), "compared to ", mean_target[component], "\n")
}

