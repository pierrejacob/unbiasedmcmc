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

load(file = "baseball.tuning.RData")
# first look at meeting times
meetingtime.list$nsamples
meetingtimes <- sapply(X = meetingtime.list$c_chains, FUN = function(x) x$meetingtime)
table(meetingtimes)
# So we can pick k = 4, and it seems like most of the time the meeting time would occur before k

# Is it a good choice?
tuning.k.list$nsamples
g <- qplot(x = tuning.k.list$ks, y = 1/(tuning.k.list$cost * tuning.k.list$v), geom = "line") + scale_y_log10()
g <- g + xlab("k") + ylab("efficiency")
g

# then we pick = 10 * k... is it a good choice?
g <- qplot(x = tuning.m.list$Ks, y = 1/(tuning.m.list$cost * tuning.m.list$v), geom = "line") + geom_point()
g <- g + xlab("m") + ylab("efficiency")
g

# OK so that seems to work nicely.
# Now we want to compare the variance with that of perfect samplers
# And that of MCMC

load(file = "baseball.mcmc.RData")
library(coda)
mcmc_var <- spectrum0(chain[1000:niterations,3])$spec
# MCMC asymptotic variance for the estimator of the mean:
mcmc_var
# For the variance of a perfect sampler, we use our coupled chains to estimate the posterior variance
load(file = "baseball.c_chain.RData")
mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_[[irep]], k = k, K = K)
}
square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_[[irep]], h = function(x) x^2, k = k, K = K)
}
estimators <- sapply(mean_estimators, function(x) x[3])
s_estimators <- sapply(square_estimators, function(x) x[3])
posterior_var <- mean(s_estimators) - mean(estimators)^2

## variance going to the variance of perfect sampler
g <- qplot(x = tuning.k.list$ks, y = tuning.k.list$v, geom = "line") + geom_point()
g <- g + scale_y_log10(breaks = c(round(posterior_var, 6), 1e-4, 1e-3, 1e-2), limits = c(1e-5, 1e-2))
g <- g + geom_hline(yintercept = posterior_var, linetype = 2) + xlab("k") + ylab("variance")
g
ggsave(filename = "baseball.tuningk.variance.pdf", plot = g, width = 5, height = 5)

# efficiency
g <- qplot(x = tuning.k.list$ks, y = 1/(tuning.k.list$cost * tuning.k.list$v), geom = "line") + scale_y_log10()
g <- g + xlab("k") + ylab("efficiency")  + geom_point()
g
ggsave(filename = "baseball.tuningk.efficiency.pdf", plot = g, width = 5, height = 5)

## variance of average estimator going to that of MCMC estimator
g <- qplot(x = tuning.m.list$Ks, y = tuning.m.list$v, geom = "line") + geom_point()
g <- g + scale_y_continuous(breaks = c(1e-6, 1e-5, round(posterior_var, 6)))
g <- g + xlab("m") + ylab("variance") + scale_x_continuous(breaks = c(4, 40, 80, 150, 200))
g <- g + geom_line(aes(y = mcmc_var / (tuning.m.list$Ks-tuning.m.list$k+1)), colour = "red", linetype = 2)
g
ggsave(filename = "baseball.tuningm.variance.pdf", plot = g, width = 5, height = 5)

load("baseball.c_chain.RData")
cost <- mean(sapply(c_chains_, function(x) x$iteration))
estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_[[irep]], h = function(x) x[3], k = k, K = K)
}
variance <- var(unlist(estimators))
efficiency_ <- 1/(cost * variance)
g <- qplot(x = tuning.m.list$Ks, y = 1/(tuning.m.list$v*tuning.m.list$cost), geom = "line") + geom_point()
g <- g + xlab("m") + ylab("efficiency")
g <- g + scale_x_continuous(breaks = c(4, 40, 80, 150, 200))
g <- g + geom_hline(yintercept = efficiency_, linetype = 2)
g
ggsave(filename = "baseball.tuningm.efficiency.pdf", plot = g, width = 5, height = 5)

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
table(meetingtime)

# g <- qplot(x = meetingtime, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time")
# g <- g + scale_x_continuous(breaks = c(2,3,4))
# g

x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g <- g + scale_x_continuous(breaks = c(2,3,4))
g
# ggsave(filename = "baseball.meetingtime.pdf", plot = g, width = 7, height = 7)

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

nclass <- 40
hist_breaks <- seq(0.37, 0.415, length.out = nclass)
hist_mcmc <- hist(chain[1000:niterations,3], breaks = hist_breaks, plot = F)
hist_mcmc$density
hist_mcmc$mids
# histogram
histogram1 <- histogram_c_chains(c_chains_, 3, k, K, breaks = hist_breaks)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[1])) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red", alpha = 0.5)
g1
ggsave(filename = "baseball.histogram.pdf", plot = g1, width = 5, height = 5)
#

hist_breaks <- seq(0.1, 0.4, length.out = nclass)
ch1 <- chain[1000:niterations,1]
hist_mcmc <- hist(ch1[ch1>0.1 & ch1 < 0.4], breaks = hist_breaks, plot = F)
# hist_mcmc <- hist(chain[1000:niterations,1], nclass = nclass, plot = F)
hist_mcmc$density
hist_mcmc$mids
histogram1 <- histogram_c_chains(c_chains_, 1, k, K, breaks = hist_mcmc$breaks)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(mu)) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red", alpha = 0.5)
g1
ggsave(filename = "baseball.histogram.mu.pdf", plot = g1, width = 5, height = 5)
