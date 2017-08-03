# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
setwd("~/Dropbox/PolyaGammaResults/reproduce/")
#
# Plummer example
##
dimension <- 2

filename <- "plummer.results.RData"
load(file = filename)
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)

# ##
K <- 1000
load(file = filename)
nsamples <- length(c_chains_)
#
k <- 500

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], k = k, K = K)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], h = function(x) x^2, k = k, K = K)
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
sum(sapply(c_chains_2, function(x) x$iteration))
hist(meetingtime)

x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g

ggsave(filename = "plummer.meetingtimes.pdf", plot = g, width = 7, height = 7)

k <- 200
K <- 5*k

load(filename)
sum(sapply(c_chains_continued_2, function(x) x$iteration))


mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_2[[irep]], k = k, K = K)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_2[[irep]], h = function(x) x^2, k = k, K = K)
}

cross_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_2[[irep]], h = function(x) x[1] * x[2], k = k, K = K)
}

est_mean <- rep(0, dimension)
est_var <- rep(0, dimension)

for (component in 1:dimension){
  estimators <- sapply(mean_estimators, function(x) x[component])
  cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
  s_estimators <- sapply(square_estimators, function(x) x[component])
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
  cat("estimated variance: ", mean(s_estimators) - mean(estimators)^2, "\n")
  est_mean[component] <- mean(estimators)
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
}
c_estimators <- sapply(cross_estimators, function(x) x[1])
cat("estimated covariance: ", mean(c_estimators) - prod(est_mean), "\n")

# histogram
histogram1 <- histogram_c_chains(c_chains_continued_2, 1, k, K, nclass = 50)
g1 <- plot_histogram(histogram1, with_bar = F) + xlab(expression(theta[2.1])) + ylab("density")
g1
ggsave("plummer.histogram1.pdf", plot = g1, width = 7, height = 7)

histogram2 <- histogram_c_chains(c_chains_continued_2, 2, k, K, nclass = 50)
g2 <- plot_histogram(histogram2, with_bar = F) + xlab(expression(theta[2.2])) + ylab("density")
g2
ggsave("plummer.histogram2.pdf", plot = g2, width = 7, height = 7)

### exact cut distribution from tedious parallel MCMC
load(file = "plummer.mcmc.RData")

hist_mcmc <- hist(theta2s[,1], breaks = histogram1$breaks, plot = F)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[2.1])) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
g1
ggsave("plummer.histogram1.pdf", plot = g1, width = 7, height = 7)

hist_mcmc <- hist(theta2s[,2], breaks = histogram2$breaks, plot = F)
g2 <- plot_histogram(histogram2, with_bar = T) + xlab(expression(theta[2.2])) + ylab("density")
g2 <- g2 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
g2
ggsave("plummer.histogram2.pdf", plot = g2, width = 7, height = 7)
