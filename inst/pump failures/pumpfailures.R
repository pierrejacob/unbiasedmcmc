# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
setwd("~/Dropbox/PolyaGammaResults/pumpfailures/")

#

##  This example is about failures of nuclear pumps. It's classic (e.g. Example 10.17 in Robert & Casella Monte Carlo Statistical Methods)
## It's used as a (laborious) example in Murdoch and Green's perfect samplers paper and also Reutter and Johnson 1995
## about using coupled chains to monitor MCMC convergence

# The data:
# number of failures
s <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
# times
t <- c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
# the model says s_k ~ Poisson(lambda_k * t_k), for k = 1,...,10
ndata <- 10
# and lambda_k ~ Gamma(alpha,beta), beta ~ Gamma(gamma, delta)
alpha <- 1.802
gamma <- 0.01
delta <- 1
# full conditionasl:
# lambda_k given rest: Gamma(alpha + s_k, beta + t_k)
# beta given rest: Gamma(gamma + 10*alpha, delta + sum_{k=1}^10 lambda_k)

single_kernel <- function(current_state){
  lambda <- current_state[1:ndata]
  beta <- current_state[ndata+1]
  for (k in 1:ndata){
    lambda[k] <- rgamma(1, shape = alpha + s[k], rate = beta + t[k])
  }
  beta <- rgamma(1, shape = gamma + 10*alpha, rate = delta + sum(lambda))
  return(c(lambda, beta))
}

rinit <- function(){
  return(rep(1, ndata+1))
}
niterations <- 5e5
# chain <- matrix(nrow = niterations, ncol = ndata+1)
# chain[1,] <- rinit()
# for (iteration in 2:niterations){
#   chain[iteration,] <- single_kernel(chain[iteration-1,])
# }
# save(chain, file = "pump.mcmc.RData")
load("pump.mcmc.RData")
hist(chain[,ndata+1])



## coupling of Gamma distributions
rgamma_coupled <- function(alpha, beta1, beta2){
  if (beta1 < beta2){
    # use Gamma(alpha, beta1) as proposal
    Z <- rgamma(1, alpha, beta1)
    logratio <- -(beta2 - beta1) * Z
    if (log(runif(1)) < logratio){
      return(list(Z1 = Z, Z2 = Z))
    } else {
      return(list(Z1 = Z, Z2 = rgamma(1, alpha, beta2)))
    }
  } else {
    Z <- rgamma(1, alpha, beta2)
    logratio <- -(beta1 - beta2) * Z
    if (log(runif(1)) < logratio){
      return(list(Z1 = Z, Z2 = Z))
    } else {
      return(list(Z1 = rgamma(1, alpha, beta1), Z2 = Z))
    }
  }
}

# alpha <- 7.5
# beta1 <- 2.3
# beta2 <- 1.1
# nsamples <- 1e5
# samp <- matrix(nrow = nsamples, ncol = 2)
# for(i in 1:nsamples){
#   x <- rgamma_coupled(alpha, beta1, beta2)
#   samp[i,] <- c(x$Z1, x$Z2)
# }
# hist(samp[,1], prob = T, nclass = 300)
# curve(dgamma(x, alpha, beta1), add=T, col = "red", lwd = 2)
# hist(samp[,2], prob = T, nclass = 300)
# curve(dgamma(x, alpha, beta2), add=T, col = "red", lwd = 2)


coupled_kernel <- function(current_state1, current_state2, ...){
  lambda1 <- current_state1[1:ndata]
  beta1 <- current_state1[ndata+1]
  lambda2 <- current_state2[1:ndata]
  beta2 <- current_state2[ndata+1]
  for (k in 1:ndata){
    x <- rgamma_coupled(alpha = alpha + s[k], beta1 = beta1 + t[k], beta2 = beta2 + t[k])
    lambda1[k] <- x$Z1
    lambda2[k] <- x$Z2
  }
  x <- rgamma_coupled(alpha =  gamma + 10*alpha, beta1 = delta + sum(lambda1), beta2 = delta + sum(lambda2))
  beta1 <- x$Z1
  beta2 <- x$Z2
  return(list(chain_state1 = c(lambda1, beta1), chain_state2 = c(lambda2, beta2)))
}


nsamples <- 10000
# c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#   coupled_chains(single_kernel, coupled_kernel, rinit)
# }
# save(c_chains_, file = "pump.cchains.RData")
load("pump.cchains.RData")

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)

x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g
ggsave(filename = "pump.meetingtime.pdf", plot = g, width = 7, height = 7)


nsamples <- 10000
k <- 20
K <- 100
# c_chains_2 <-  foreach(irep = 1:nsamples) %dorng% {
#   coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
# }
# save(c_chains_2, file = "pump.cchains2.RData")
load("pump.cchains2.RData")

meetingtime <- sapply(c_chains_2, function(x) x$meetingtime)
summary(meetingtime)

iterations <- sapply(c_chains_2, function(x) x$iteration)
sum(iterations)

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], k = k, K = K)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], h = function(x) x^2, k = k, K = K)
}

est_mean <- rep(0, ndata+1)
est_var <- rep(0, ndata+1)
for (component in 1:(ndata+1)){
  cat("component ", component, "\n")
  estimators <- sapply(mean_estimators, function(x) x[component])
  cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
  s_estimators <- sapply(square_estimators, function(x) x[component])
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
  cat("estimated variance: ", mean(s_estimators) - mean(estimators)^2, "\n")
  est_mean[component] <- mean(estimators)
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
}

# histogram
histogram1 <- histogram_c_chains(c_chains_2, ndata+1, k, K, nclass = 100)
# g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[1])) + ylab("density")
# g1
# # looks similar to
# qplot(x = chain[,ndata+1], geom = "blank") + geom_histogram(aes(y = ..density..))

hist_mcmc <- hist(chain[1000:niterations,ndata+1], breaks = histogram1$breaks, plot = F)
# histogram
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[1])) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
g1
ggsave(filename = "pump.histogram.pdf", plot = g1, width = 7, height = 7)


