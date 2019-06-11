# load packages
library(unbiasedmcmc)
library(latex2exp)
library(dplyr)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())

#
# Plummer's example
dimension <- 2

filename <- "plummer.tuning.RData"
load(file = filename)
meetingtime <- sapply(c_chains_1, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)
# from which we have found
k
m

# ##
load(file = "plummer.results.RData")
nsamples
k
m

meetingtime <- sapply(c_chains_2, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)
#

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], k = k, m = m)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_2[[irep]], h = function(x) x^2, k = k, m = m)
}

est_mean <- rep(0, dimension)
est_var <- rep(0, dimension)
for (component in 1:dimension){
  estimators <- sapply(mean_estimators, function(x) x[component])
  est_mean[component] <- mean(estimators)
  cat("estimated mean: ", est_mean[component], "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
  s_estimators <- sapply(square_estimators, function(x) x[component])
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
  cat("estimated variance: ", est_var[component], "\n")
}

nsamples <- length(c_chains_2)
meetingtime <- sapply(c_chains_2, function(x) x$meetingtime)
summary(meetingtime)

# ggsave(filename = "plummer.meetingtimes.pdf", plot = g, width = 7, height = 7)
gmeetingtime <- qplot(x = meetingtime, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time")
gmeetingtime
ggsave(filename = "plummer.meetingtimes.pdf", plot = gmeetingtime, width = 5, height = 5)

sum(sapply(c_chains_2, function(x) x$iteration)) / nsamples

### cut distribution from tedious parallel MCMC
histogram1 <- histogram_c_chains(c_chains_2, 1, k, m, nclass = 35)
histogram2 <- histogram_c_chains(c_chains_2, 2, k, m, nclass = 30)
load(file = "plummer.mcmc.RData")

hist_mcmc <- hist(theta2s[,1], breaks = histogram1$breaks, plot = F)
# hist_mcmc <- hist(theta2s[,1], plot = F)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(TeX("$\\theta_{2,1}$")) + ylab("density")
g1 <- g1 + geom_line(data=data.frame(x = hist_mcmc$mids, y = hist_mcmc$density), aes(x = x, y = y, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL), colour = "red")
g1 <- g1 + scale_x_continuous(breaks = c(-2.5, -2, -1.5))
g1
ggsave("plummer.histogram1.pdf", plot = g1, width = 5, height = 5)

hist_mcmc <- hist(theta2s[,2], breaks = histogram2$breaks, plot = F)
g2 <- plot_histogram(histogram2, with_bar = T) + xlab(TeX("$\\theta_{2,2}$")) + ylab("density")
g2 <- g2 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
g2 <- g2 + scale_x_continuous(breaks = c(10, 15, 20, 25))
g2
ggsave("plummer.histogram2.pdf", plot = g2, width = 5, height = 5)
