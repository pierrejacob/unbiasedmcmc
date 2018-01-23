
# Load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores()-2)
setwd("~/Dropbox/UnbiasedMCMCResults/germancredit/")
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
logistic_setting <- logistic_precomputation(Y, X, b, B)


# MCMC transition kernels
single_kernel <- function(chain_state, logistic_setting){
  zs <- abs(xbeta(logistic_setting$X, t(chain_state)))
  w <- rpg(logistic_setting$n, h=1, z=zs)
  res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(chain_state)
}

coupled_kernel <- function(chain_state1, chain_state2, logistic_setting, return_ws=FALSE){
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  if(!return_ws){
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2)))
  } else {
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2),w1=ws$w1,w2=ws$w2))
  }
}

rinit <- function(){
  t(fast_rmvnorm(1, mean = b, covariance = B))
}


# Behavior from a single run
current_value1 <- rinit()
current_value2 <- rinit()
niterations <- 1000

chain1_w <- chain2_w <- matrix(ncol=n, nrow=niterations)
chain1_beta <- chain2_beta <- matrix(ncol=p, nrow=niterations)
chain1_beta[1,] <- current_value1
chain2_beta[1,] <- current_value2

for(t in 2:niterations){
  current_value <- coupled_kernel(current_value1, current_value2, logistic_setting, return_ws=TRUE)
  chain1_w[t-1,] <- current_value$w1
  chain2_w[t-1,] <- current_value$w2
  chain1_beta[t,] <- current_value1 <- current_value$chain_state1
  chain2_beta[t,] <- current_value2 <- current_value$chain_state2
  if(all(current_value1==current_value2)) break
}
meetingtime=t

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
nsamples <- 1000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit, logistic_setting = logistic_setting)
}

meetingtime.list <- list(c_chains = c_chains_, nsamples = nsamples)

#meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
#table(meetingtime)
#quantile(meetingtime, .95)

save(df_beta, df_w, meetingtime.list, file = "germancredit.tuning.RData")
load(file = "germancredit.tuning.RData")

# Histogram of variable estimates
k <- 110
m <- 1100

nsamples <- 1000

c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit, m=m, logistic_setting=logistic_setting)
}
save(nsamples, k, m, c_chains_, file = "germancredit.c_chain.RData")


# Efficiency for a grid of k values
idx1 <- which('Instalment.per.cent' == colnames(X))
idx2 <- which('Duration.in.Current.address' == colnames(X))
idx <- idx2
colnames(X)[idx]

nsamples <- 1000
ks <- seq(0, 250, by=5)
cost <- rep(0, length(ks))
v <- rep(0, length(ks))
for (ik in 1:length(ks)){
  k <- ks[ik]
  cat("k = ", k, "\n")
  m <- k
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(single_kernel, coupled_kernel, rinit, m = m, logistic_setting=logistic_setting)
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











