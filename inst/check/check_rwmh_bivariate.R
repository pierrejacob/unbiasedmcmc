library(unbiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

## test MH kernels on a univariate target
## define target distribution
dimension <- 2
Sigma_pi <- diag(1, dimension, dimension)
target <- function(x){
  evals <- rep(0, 2)
  evals[1] <- log(0.5) + fast_dmvnorm(matrix(x, nrow = 1), mean = c(0, 0), covariance = Sigma_pi)
  evals[2] <- log(0.5) + fast_dmvnorm(matrix(x, nrow = 1), mean = c(4, 0), covariance = Sigma_pi)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

Sigma_proposal <- diag(2, dimension, dimension)
Sigma_proposal_chol <- chol(Sigma_proposal)
Sigma_proposal_chol_inv <- solve(chol(Sigma_proposal))

kernels <- get_mh_kernels(target, Sigma_proposal)
single_kernel <- kernels$single_kernel
coupled_kernel <- kernels$coupled_kernel
#
rinit <- function(){
  x <- fast_rmvnorm(1, rep(5, dimension), diag(5, dimension, dimension))
  return(list(chain_state = x, current_pdf = target(x)))
}

### Test of the single kernel
niterations <- 50000
state <- rinit()
chain <- matrix(ncol=dimension, nrow=niterations)
chain[1,] <- state$chain_state
for (t in 2:niterations){
  state <- single_kernel(state)
  chain[t,] <- state$chain_state
}
it <- floor(seq(from = 1, to  = niterations, length.out = 1000))
qplot(x=it, y=chain[it,1], geom = "line") + ylab("X") + xlab("iteration")
qplot(x=it, y=chain[it,2], geom = "line") + ylab("X") + xlab("iteration")

hist(chain[1000:niterations,1], nclass = 100, prob = TRUE)
## numerical integration to obtain marginal
marg1 <- function(x1){
  return(integrate(f = function(x) sapply(x, function(v) exp(target(matrix(c(x1,v), nrow = 1)))), lower = -10, upper = 10)$value)
}
curve(sapply(x, function(y) marg1(y)), add = T)

hist(chain[1000:niterations,2], nclass = 100, prob = TRUE)
marg2 <- function(x2){
  return(integrate(f = function(x) sapply(x, function(v) exp(target(matrix(c(v,x2), nrow = 1)))), lower = -10, upper = 10)$value)
}
curve(sapply(x, function(y) marg2(y)), add = T)

# Markov kernel of the coupled chain
niterations <- 50000
state1 <- rinit()
state2 <- rinit()
chain1 <- matrix(ncol=dimension, nrow=niterations)
chain2 <- matrix(ncol=dimension, nrow=niterations)
chain1[1,] <- state1$chain_state
chain2[1,] <- state2$chain_state
meetingtime <- Inf
for (t in 2:niterations){
  coupled_result <- coupled_kernel(state1, state2)
  state1 <- coupled_result$state1
  state2 <- coupled_result$state2
  if (coupled_result$identical && is.infinite(meetingtime)){
    meetingtime <- t
  }
  chain1[t,] <- state1$chain_state
  chain2[t,] <- state2$chain_state
}


hist(chain1[1000:niterations,1], nclass = 100, prob = TRUE)
hist(chain2[1000:niterations,1], nclass = 100, prob = TRUE, add = TRUE, col = rgb(1,0,0,0.5))
## numerical integration to obtain marginal
marg1 <- function(x1){
  return(integrate(f = function(x) sapply(x, function(v) exp(target(matrix(c(x1,v), nrow = 1)))), lower = -10, upper = 10)$value)
}
curve(sapply(x, function(y) marg1(y)), add = T)

hist(chain1[1000:niterations,2], nclass = 100, prob = TRUE)
hist(chain2[1000:niterations,2], nclass = 100, prob = TRUE, add = TRUE, col = rgb(1,0,0,0.5))
marg2 <- function(x2){
  return(integrate(f = function(x) sapply(x, function(v) exp(target(matrix(c(v,x2), nrow = 1)))), lower = -10, upper = 10)$value)
}
curve(sapply(x, function(y) marg2(y)), add = T)



res <- sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = 3)
res$meetingtime
res$iteration
H_bar(res, h = function(x) x, k = 3, m = 3)
# #
qplot(x = 0:res$iteration, y = res$samples1[,1], geom = "line") +
  geom_line(aes(x = 0:(res$iteration-1), y = res$samples2[,1]), col = rgb(1,0,0,0.5)) +
  geom_line(aes(x = 1:res$iteration, y = res$samples2[,1]), col = rgb(1,0,0,1), linetype = 2)

### distribution of meeting times
nsamples <- 500
meetingtime <- foreach(irep = 1:nsamples, .combine = c) %dorng% {
  c_chains <- sample_meetingtime(single_kernel, coupled_kernel, rinit)
  c_chains$meetingtime
}
summary(meetingtime)
hist(meetingtime)


