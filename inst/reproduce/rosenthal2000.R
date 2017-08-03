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

niterations <- 5e5
# chain <- matrix(nrow = niterations, ncol = ndata+2)
# chain[1,] <- rinit()
# for (iteration in 2:niterations){
#   chain[iteration,] <- single_kernel(chain[iteration-1,])
# }
# save(chain, file = "rosenthal2000.chain.RData")
load(file = "rosenthal2000.chain.RData")

nsamples <- 10000
# c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#   coupled_chains(single_kernel, coupled_kernel, rinit)
# }
# save(c_chains_, file = "rosenthal2000.c_chain.RData")
load(file = "rosenthal2000.c_chain.RData")

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
tabulate(meetingtime)
# g <- qplot(x = meetingtime, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time")
# g <- g + scale_x_continuous(breaks = c(2,3,4))
# g
x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g <- g + scale_x_continuous(breaks = c(2,3,4))
g
ggsave(filename = "rosenthal2000.meetingtime.pdf", plot = g, width = 7, height = 7)

nsamples <- 10000
k <- 4
K <- 50
# c_chains_2 <-  foreach(irep = 1:nsamples) %dorng% {
#   coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
# }
# save(c_chains_2, file = "rosenthal2000.c_chain2.RData")
load(file = "rosenthal2000.c_chain2.RData")

meetingtime <- sapply(c_chains_2, function(x) x$meetingtime)
tabulate(meetingtime)

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], k = k, K = K)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], h = function(x) x^2, k = k, K = K)
}

# cross_estimators <-  foreach(irep = 1:nsamples) %dorng% {
#   H_bar(c_chains_2[[irep]], h = function(x) x[1] * x[2], k = k, K = K)
# }

est_mean <- rep(0, ndata+2)
est_var <- rep(0, ndata+2)
for (component in 1:(ndata+2)){
  cat("component ", component, "\n")
  estimators <- sapply(mean_estimators, function(x) x[component])
  cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
  s_estimators <- sapply(square_estimators, function(x) x[component])
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
  cat("estimated variance: ", mean(s_estimators) - mean(estimators)^2, "\n")
  est_mean[component] <- mean(estimators)
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
}

hist_breaks <- seq(0.37, 0.415, length.out = 100)
hist_mcmc <- hist(chain[1000:niterations,3], breaks = hist_breaks, plot = F)
hist_mcmc$density
hist_mcmc$mids
# histogram
histogram1 <- histogram_c_chains(c_chains_2, 3, k, K, breaks = hist_breaks)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[1])) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
g1
ggsave(filename = "rosenthal2000.histogram.pdf", plot = g1, width = 7, height = 7)
# looks similar to
# qplot(x = chain[,3], geom = "blank") + geom_histogram(aes(y = ..density..)) + xlim(0.36, 0.425)


