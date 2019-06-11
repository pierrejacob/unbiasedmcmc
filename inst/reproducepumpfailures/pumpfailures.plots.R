
# load packages
library(unbiasedmcmc)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

ndata <- 10

load("pump.tuning.RData")
# first look at meeting times
nsamples <- meetingtime.list$nsamples
meetingtimes <- sapply(X = meetingtime.list$c_chains, FUN = function(x) x$meetingtime)
table(meetingtimes)

# So we can pick k = 7, and it seems like most of the time the meeting time would occur before k

# Is it a good choice?
tuning.k.list$nsamples
g <- qplot(x = tuning.k.list$ks, y = 1/(tuning.k.list$cost * tuning.k.list$v), geom = "line") + scale_y_continuous(breaks = c(0.2, 0.3, 0.4))
g <- g + xlab("k") + ylab("efficiency")
g

# then we pick = 10 * k... is it a good choice?
g <- qplot(x = tuning.m.list$ms, y = 1/(tuning.m.list$cost * tuning.m.list$v), geom = "line") + geom_point()
g <- g + xlab("m") + ylab("efficiency")
g

load("pump.mcmc.RData")
library(coda)
mcmc_var <- spectrum0(chain[1000:niterations,ndata+1])$spec
# MCMC asymptotic variance for the estimator of the mean:
mcmc_var

load("pump.c_chain.RData")
mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_[[irep]], k = k, m = m)
}
square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_[[irep]], h = function(x) x^2, k = k, m = m)
}
estimators <- sapply(mean_estimators, function(x) x[ndata+1])
s_estimators <- sapply(square_estimators, function(x) x[ndata+1])
posterior_var <- mean(s_estimators) - mean(estimators)^2

## variance going to the variance of perfect sampler
g <- qplot(x = tuning.k.list$ks, y = tuning.k.list$v, geom = "line") + geom_point()
g <- g + geom_hline(yintercept = posterior_var, linetype = 2) + xlab("k") + ylab("variance")
g
ggsave(filename = "pump.tuningk.variance.pdf", plot = g, width = 5, height = 5)

# efficiency
g <- qplot(x = tuning.k.list$ks, y = 1/(tuning.k.list$cost * tuning.k.list$v), geom = "line") + geom_point()
g <- g + xlab("k") + ylab("efficiency")
g
ggsave(filename = "pump.tuningk.efficiency.pdf", plot = g, width = 5, height = 5)

## variance of average estimator going to that of MCMC estimator
g <- qplot(x = tuning.m.list$ms, y = tuning.m.list$v, geom = "line") + geom_point()
g <- g + xlab("m") + ylab("variance") + scale_x_continuous(breaks = c(4, 20, 40, 80, 150, 200))
g <- g + geom_line(aes(y = mcmc_var / (tuning.m.list$ms-tuning.m.list$k+1)), colour = "red", linetype = 2)
g
ggsave(filename = "pump.tuningm.variance.pdf", plot = g, width = 5, height = 5)

# efficiency of our choice k = 7, K = 70
load("pump.c_chain.RData")
cost <- mean(sapply(c_chains_, function(x) x$iteration))
estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_[[irep]], h = function(x) x[ndata+1], k = k, m = m)
}
variance <- var(unlist(estimators))
efficiency_ <- 1/(cost * variance)

g <- qplot(x = tuning.m.list$ms, y = 1/(tuning.m.list$v*tuning.m.list$cost), geom = "line") + geom_point()
g <- g + xlab("m") + ylab("efficiency")
g <- g + scale_x_continuous(breaks = c(4, 20, 40, 80, 150, 200))
g <- g + geom_hline(yintercept = efficiency_, linetype = 2)
g
ggsave(filename = "pump.tuningm.efficiency.pdf", plot = g, width = 5, height = 5)

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
table(meetingtime)


meetingtime <- sapply(meetingtime.list$c_chains, function(x) x$meetingtime)
x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g <- g + scale_x_continuous(breaks = 2:10)
g
g <- qplot(x = meetingtime, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time")
g <- g + scale_x_continuous(breaks = 2:10)
g
ggsave(filename = "pump.meetingtime.pdf", plot = g, width = 5, height = 5)
# ggsave(filename = "pump.meetingtime.pdf", plot = g, width = 7, height = 7)
k
m

# mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
#   H_bar(c_chains_[[irep]], k = k, K = K)
# }
#
# square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
#   H_bar(c_chains_[[irep]], h = function(x) x^2, k = k, K = K)
# }
# est_mean <- rep(0, ndata+1)
# est_var <- rep(0, ndata+1)
# for (component in 1:(ndata+1)){
#   cat("component ", component, "\n")
#   estimators <- sapply(mean_estimators, function(x) x[component])
#   cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
#   s_estimators <- sapply(square_estimators, function(x) x[component])
#   cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
#   cat("estimated variance: ", mean(s_estimators) - mean(estimators)^2, "\n")
#   est_mean[component] <- mean(estimators)
#   est_var[component] <- mean(s_estimators) - est_mean[component]^2
# }

# histogram
nclass <- 27
histogram1 <- histogram_c_chains(c_chains_, ndata+1, k, m, nclass = nclass)
niterations <- nrow(chain)
hist_mcmc <- hist(chain[1000:niterations,ndata+1], nclass = 100, plot = FALSE)

g1 <- plot_histogram(histogram1, with_bar = TRUE) + xlab(expression(beta)) + ylab("density")
g1 <- g1 + geom_line(data = data.frame(x = hist_mcmc$mids, y = hist_mcmc$density), aes(x = x, y = y, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL), colour = "red")
g1
ggsave(filename = "pump.histogram.pdf", plot = g1, width = 5, height = 5)
