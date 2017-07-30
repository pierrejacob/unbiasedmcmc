# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
setwd("~/Dropbox/PolyaGammaResults/rosenthal2000/")
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
  As <- rigamma_coupled(alpha = a + 0.5 * (ndata-1), beta1 = b + 0.5 * sum((theta1 - theta_bar1)^2),
                        beta2 = b + 0.5 * sum((theta2 - theta_bar2)^2))
  A1 <- As$Z1
  A2 <- As$Z2
  # update of mu given rest
  mus <- gaussian_max_coupling(theta_bar1, theta_bar2, matrix(A1/ndata, 1, 1), matrix(A2/ndata, 1, 1))
  mu1 <- mus[,1]
  mu2 <- mus[,2]
  # # update of each theta_i
  thetas <- gaussian_max_coupling((mu1 * V + Y * A1) / (V + A1),
                                  (mu2 * V + Y * A2) / (V + A2),
                                  diag(A1 * V / (V + A1), ndata, ndata),
                                  diag(A2 * V / (V + A2), ndata, ndata))
  theta1 <- thetas[,1]
  theta2 <- thetas[,2]
  return(list(chain_state1 = c(mu1, A1, theta1), chain_state2 = c(mu2, A2, theta2)))
}

nsamples <- 10000
ks <- 2:10
# cost <- rep(0, length(ks))
# v <- rep(0, length(ks))
# for (ik in 1:length(ks)){
#   k <- ks[ik]
#   cat("k = ", k, "\n")
#   K <- k
#   c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#     coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
#   }
#   estimators <-  foreach(irep = 1:nsamples) %dorng% {
#     H_bar(c_chains_[[irep]], h = function(x) x[3], k = k, K = K)
#   }
#   cost[ik] <- sum(sapply(c_chains_, function(x) x$iteration)) / nsamples
#   v[ik] <- var(unlist(estimators))
# }
# save(cost, v, file = "rosenthal2000.various_k.RData")
load(file = "rosenthal2000.various_k.RData")

# estimate posterior variance of component 3
load(file = "rosenthal2000.c_chain2.RData")
mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], k = 4, K = 50)
}
square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], h = function(x) x^2, k = 4, K = 50)
}
estimators <- sapply(mean_estimators, function(x) x[3])
s_estimators <- sapply(square_estimators, function(x) x[3])
post_var <- mean(s_estimators) - mean(estimators)^2

g <- qplot(x = ks, y = v, geom = "line") + scale_x_continuous(breaks = ks) + scale_y_log10(breaks = c(1e-5, round(post_var, 6), 1e-4, 1e-3), limits = c(1e-5, 1e-2))
g <- g + geom_hline(yintercept = post_var, linetype = 2) + xlab("k") + ylab("variance")
g
ggsave(filename = "rosenthal2000.variance.k.pdf", plot = g, width = 7, height = 7)

# qplot(x = ks, y = cost, geom = "line") + scale_x_continuous(breaks = ks)

k <- 4
Ks <- floor(seq(from = 4, to = 50, length.out = 10))
# cost <- rep(0, length(Ks))
# v <- rep(0, length(Ks))
# for (iK in 1:length(Ks)){
#   K <- Ks[iK]
#   cat("K = ", K, "\n")
#   c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#     coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
#   }
#   estimators <-  foreach(irep = 1:nsamples) %dorng% {
#     H_bar(c_chains_[[irep]], h = function(x) x[3], k = k, K = K)
#   }
#   cost[iK] <- sum(sapply(c_chains_, function(x) x$iteration)) / nsamples
#   v[iK] <- var(unlist(estimators))
# }
# save(cost, v, file = "rosenthal2000.various_ktoK.RData")
load(file = "rosenthal2000.various_ktoK.RData")

g <- qplot(x = Ks, y = v, geom = "line") + scale_x_continuous(breaks = Ks) + scale_y_continuous(breaks = c(1e-6, 1e-5, round(post_var, 6)))
g <- g + xlab("m") + ylab("variance")
g <- g + geom_line(aes(y = mcmc_var$spec / (Ks-k+1)), colour = "red", linetype = 2)
g
ggsave(filename = "rosenthal2000.variance.ktoK.pdf", plot = g, width = 7, height = 7)

library(mcmcse)
load("rosenthal2000.chain.RData")
niterations <- nrow(chain)
mcmcse::mcse(chain[1000:niterations,3])$se
?mcmcse::mcse
library(mcmc)
library(MCMCpack)
mcmc_var <- spectrum0(chain[1000:niterations,3])

