
# Load packages
library(debiasedmcmc)
library(dplyr)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
#
#  This example is about the Polya Gamma Gibbs sampler for logistic regression models, as applied to the German credit data of Lichman 2013.
#

# Data
# data <- read.csv('german_credit.csv')
# Y <- data[,"Creditability"]
# x.categorical <- c('Account.Balance', 'Payment.Status.of.Previous.Credit', 'Purpose', 'Value.Savings.Stocks',
#                    'Length.of.current.employment', 'Sex...Marital.Status', 'Guarantors', 'Most.valuable.available.asset',
#                    'Concurrent.Credits', 'Type.of.apartment', 'Occupation', 'Telephone', 'Foreign.Worker')
# x.quant <- c('Duration.of.Credit..month.', 'Credit.Amount', 'Instalment.per.cent', 'Duration.in.Current.address',
#              'Age..years.', 'No.of.Credits.at.this.Bank', 'No.of.dependents')
# for(x in x.categorical){
#   data[,x] = as.factor(data[,x])
# }
# fmla <- paste('~',paste(c(x.quant,x.categorical),collapse ='+'))
# X <- model.matrix(formula(fmla), data=data)
data(germancredit)
n <- nrow(X)
p <- ncol(X)


# Prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
logistic_setting <- logisticregression_precomputation(Y, X, b, B)


# MCMC transition kernels
# single_kernel <- function(chain_state, logistic_setting){
#   zs <- abs(xbeta(logistic_setting$X, t(chain_state)))
#   w <- rpg(logistic_setting$n, h=1, z=zs)
#   res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
#   chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
#   return(chain_state)
# }

single_kernel <- function(state){
  beta <- state$chain_state[1:p]
  zs <- abs(logisticregression_xbeta(logistic_setting$X, beta))
  w <- BayesLogit::rpg(logistic_setting$n, h=1, z=zs)
  res <- logisticregression_m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- c(fast_rmvnorm_chol(1, res$m, res$Cholesky), w)
  return(list(chain_state = chain_state))
}

coupled_kernel <- function(state1, state2){
  beta1 <- state1$chain_state[1:p]
  beta2 <- state2$chain_state[1:p]
  ws <- sample_w(beta1, beta2, logistic_setting$X)
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  identical = FALSE
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
    identical = TRUE
  }
  return(list(state1 = list(chain_state=c(betas$beta1, ws$w1)),
              state2 = list(chain_state=c(betas$beta2, ws$w2)),
              identical = identical))
}

rinit <- function(){
  return(list(chain_state = c(fast_rmvnorm(1, mean = b, covariance = B)[1,], rep(0, n))))
}


# Behavior from a single run
state1 <- rinit()
state2 <- rinit()
niterations <- 1000
coupled_state <- list(state1 = state1, state2 = state2, identical = FALSE)
chain1_w <- chain2_w <- matrix(ncol=n, nrow=niterations)
chain1_beta <- chain2_beta <- matrix(ncol=p, nrow=niterations)
chain1_beta[1,] <- state1$chain_state[1:p]
chain2_beta[1,] <- state2$chain_state[1:p]

for(t in 2:niterations){
  coupled_state <- coupled_kernel(coupled_state$state1, coupled_state$state2)
  chain1_w[t-1,] <- coupled_state$state1$chain_state[(p+1):(p+n)]
  chain2_w[t-1,] <- coupled_state$state2$chain_state[(p+1):(p+n)]
  chain1_beta[t,] <- coupled_state$state1$chain_state[1:p]
  chain2_beta[t,] <- coupled_state$state2$chain_state[1:p]
  if(coupled_state$identical) break
}
meetingtime=t
meetingtime

# Compute distances and number met for plotting
nmet_w <- sapply(1:meetingtime, function(t) sum(chain1_w[t,]==chain2_w[t,]))
nmet_w[meetingtime] <- n
distances_w <- sapply(1:meetingtime, function(t) sqrt(sum((chain1_w[t,] - chain2_w[t,])^2)))
distances_w[meetingtime] <- 0
distances_beta <- sapply(1:meetingtime, function(t) sqrt(sum((chain1_beta[t,] - chain2_beta[t,])^2)))
df_beta <- data.frame(iter=1:meetingtime, dist_beta = distances_beta[1:meetingtime])
df_w <- data.frame(iter=1:meetingtime, dist_w = distances_w[1:meetingtime], nmet_w = nmet_w[1:meetingtime])

save(df_beta, df_w, file = "germancredit.tuning.RData")
load(file = "germancredit.tuning.RData")


# Distribution of meeting times over many runs
nsamples <- 100
meetings_1 <-  foreach(irep = 1:nsamples) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit)
}

meetingtime.list <- list(meetings_1 = meetings_1, nsamples = nsamples)

#meetingtime <- sapply(meetings_1, function(x) x$meetingtime)
#table(meetingtime)
#quantile(meetingtime, .95)

save(df_beta, df_w, meetingtime.list, file = "germancredit.tuning.RData")
load(file = "germancredit.tuning.RData")

# Histogram of variable estimates
k <- 110
m <- 1100

nsamples <- 100

c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m)
}
save(nsamples, k, m, c_chains_, file = "germancredit.c_chain.RData")


# Efficiency for a grid of k values
idx1 <- which('Instalment.per.cent' == colnames(X))
idx2 <- which('Duration.in.Current.address' == colnames(X))
idx <- idx2
colnames(X)[idx]

nsamples <- 100
ks <- seq(0, 250, by=5)
cost <- rep(0, length(ks))
v <- rep(0, length(ks))
for (ik in 1:length(ks)){
  k <- ks[ik]
  cat("k = ", k, "\n")
  m <- k
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m)
  }
  estimators <-  foreach(irep = 1:nsamples) %dorng% {
    H_bar(c_chains_[[irep]], h = function(x) x[idx], k = k, m = m)
  }
  cost[ik] <- sum(sapply(c_chains_, function(x) x$iteration)) / nsamples
  v[ik] <- var(unlist(estimators))
}

tuning.k.list <- list(nsamples = nsamples, ks = ks, cost = cost, v = v)
save(df_beta, df_w, meetingtime.list, tuning.k.list, file = "germancredit.tuning.RData")

#qplot(x = ks, y = 1 / (cost * v), geom = "line")
#k <- ks[which.max(1/(v*cost))]











