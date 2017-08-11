# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
# this example is taken from Rosenthal "Parallel Computing and Monte Carlo algorithms", 2000.
# The data are baseball player's batting averages, as in Efron and Morris 1975 or Morris 1983
# In particular, the data are listed in Morris 1983, and available as part of Rosenthal's online code:
# http://probability.ca/jeff/comp/james.c

Y <- c(0.395, 0.375, 0.355, 0.334, 0.313, 0.313, 0.291, 0.269, 0.247, 0.247, 0.224, 0.224,
       0.224, 0.224, 0.224, 0.200, 0.175, 0.148)
ndata <- length(Y)

# The data are Y_i for i in {1,...,K}, with here K = 18.
# The model specifies Y_i ~ Normal(theta_i, V), with V known and fixed to 0.00434
V <- 0.00434
# The prior is theta_i ~ Normal(mu, A), with mu ~ flat prior, and A ~ InverseGamma(a,b)
# with density proportional to exp(-b/x) x^{-a-1}; Rosenthal chooses a = -1, b = 2
a <- -1
b <- 2
# To target the posterior distribution, we follow Rosenthal and consider the following Gibbs sampler.
# A given rest: IG(a + (K-1)/2, b + (sum_i=1^K (theta_i - theta_bar)^2)/2)
# mu given rest: Normal(theta_bar, A / K)
# theta_i given rest: Normal( (mu * V + Y_i * A) / (V + A), A * V / (V + A))
# ... where theta_bar is the average of the theta_i, for i in {1,...,K}


# we store the parameters as (mu, A, theta_1, ..., theta_K) so the parameter space is of dimension 20
# the initialization comes from Rosenthal (except the initialization of mu and A which is irrelevant)
rinit <- function(){
  return(c(0,1,rep(mean(Y), ndata)))
}
# here is the code for the Gibbs sampler
single_kernel <- function(current_state, ...){
  theta <- current_state[3:(ndata+2)]
  theta_bar <- mean(current_state[3:(ndata+2)])
  # update of A given rest
  A <- rigamma(1, a + 0.5 * (ndata-1), b + 0.5 * sum((theta - theta_bar)^2))
  # update of mu given rest
  mu <- rnorm(1, theta_bar, A/ndata)
  # update of each theta_i
  theta <- rnorm(ndata, (mu * V + Y * A) / (V + A), A * V / (V + A))
  return(c(mu, A, theta))
}

# Now coupled kernel
coupled_kernel <- function(current_state1, current_state2, ...){
  theta1 <- current_state1[3:(ndata+2)]
  theta_bar1 <- mean(theta1)
  theta2 <- current_state2[3:(ndata+2)]
  theta_bar2 <- mean(theta2)
  # update of A given rest
  As <- rigamma_coupled(alpha1 = a + 0.5 * (ndata-1), alpha2 = a + 0.5 * (ndata-1),
                        beta1 = b + 0.5 * sum((theta1 - theta_bar1)^2),
                        beta2 = b + 0.5 * sum((theta2 - theta_bar2)^2))
  A1 <- As[1]
  A2 <- As[2]
  # update of mu given rest
  # mus <- gaussian_max_coupling(theta_bar1, theta_bar2, matrix(A1/ndata, 1, 1), matrix(A2/ndata, 1, 1))
  mus <- rnorm_max_coupling(theta_bar1, theta_bar2, sqrt(A1/ndata), sqrt(A2/ndata))
  mu1 <- mus[1]
  mu2 <- mus[2]
  # # update of each theta_i
  thetas <- matrix(nrow = ndata, ncol = 2)
  for (idata in 1:ndata){
    thetas[idata,] <- rnorm_max_coupling((mu1 * V + Y[idata] * A1) / (V + A1), (mu2 * V + Y[idata] * A2) / (V + A2),
                                         sqrt(A1 * V / (V + A1)), sqrt(A2 * V / (V + A2)))
  }
  theta1 <- thetas[,1]
  theta2 <- thetas[,2]
  return(list(chain_state1 = c(mu1, A1, theta1), chain_state2 = c(mu2, A2, theta2)))
}

### First, get a sample of i.i.d. meeting times
nsamples <- 1000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit)
}
meetingtime.list <- list(c_chains = c_chains_, nsamples = nsamples)
save(meetingtime.list, file = "baseball.tuning.RData")
load(file = "baseball.tuning.RData")
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
table(meetingtime)
# we have seen only times = to 2 and 3
# so we pick k = 4

#
nsamples <- 10000
ks <- 2:10
cost <- rep(0, length(ks))
v <- rep(0, length(ks))
for (ik in 1:length(ks)){
  k <- ks[ik]
  cat("k = ", k, "\n")
  K <- k
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
  }
  estimators <-  foreach(irep = 1:nsamples) %dorng% {
    H_bar(c_chains_[[irep]], h = function(x) x[3], k = k, K = K)
  }
  cost[ik] <- sum(sapply(c_chains_, function(x) x$iteration)) / nsamples
  v[ik] <- var(unlist(estimators))
}
tuning.k.list <- list(nsamples = nsamples, ks = ks, cost = cost, v = v)
save(meetingtime.list, tuning.k.list, file = "baseball.tuning.RData")

load(file = "baseball.tuning.RData")
qplot(x = ks, y = cost * v, geom = "line") + scale_y_log10()
# clearly k = 4 was a good choice !

# Then with our heuristics we would choose m = 10 * k,
# so m = 40.
# Let's try different values of m and compare.

k <- 4
Ks <- c(4, 10, 20, 40, 60, 80, 100, 150, 200)
Ks
cost <- rep(0, length(Ks))
v <- rep(0, length(Ks))
for (iK in 1:length(Ks)){
  K <- Ks[iK]
  cat("K = ", K, "\n")
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
  }
  estimators <-  foreach(irep = 1:nsamples) %dorng% {
    H_bar(c_chains_[[irep]], h = function(x) x[3], k = k, K = K)
  }
  cost[iK] <- sum(sapply(c_chains_, function(x) x$iteration)) / nsamples
  v[iK] <- var(unlist(estimators))
}
tuning.m.list <- list(nsamples = nsamples, k = k, Ks = Ks, cost = cost, v = v)
save(meetingtime.list, tuning.k.list, tuning.m.list, file = "baseball.tuning.RData")
load(file = "baseball.tuning.RData")

qplot(x = Ks, y = 1/(cost * v), geom = "line") + geom_point() + geom_vline(xintercept = 4 * 10)
# achieves nearly maximum efficiency


# library(mcmcse)
# niterations <- nrow(chain)
# mcmcse::mcse(chain[1000:niterations,3])$se
# ?mcmcse::mcse
# library(mcmc)
# library(MCMCpack)
# mcmc_var <- spectrum0(chain[1000:niterations,3])

