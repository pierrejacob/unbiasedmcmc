# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores())
#

load(file = "propensity.data.RData")
load(file = "propensity.second_data.RData")

load(file = "propensity.stage1.c_chains.RData")
meetingtime <- sapply(c_chains_continued_, function(x) x$meetingtime)
summary(meetingtime)
nsamples <- length(c_chains_continued_)
# x <- as.numeric(names(table(meetingtime)))
# y <- as.numeric(table(meetingtime)) / nsamples
# g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
# g
g <- qplot(x = meetingtime, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time")
g <- g + scale_x_continuous(breaks = 2:10)
g
ggsave(filename = "propensity.meetingtime.stage1.pdf", plot = g, width = 5, height = 5)


load(file = "propensity.stage2.c_chains.RData")
nsamples <- second_c_chains_ %>% length
meetingtime <- sapply(second_c_chains_, function(x) x$meetingtime)
summary(meetingtime)
# x <- as.numeric(names(table(meetingtime)))
# y <- as.numeric(table(meetingtime)) / nsamples
# g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
# g
g <- qplot(x = meetingtime, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time")
g <- g + scale_x_continuous(breaks = 2:10)
g
ggsave(filename = "propensity.meetingtime.stage2.pdf", plot = g, width = 5, height = 5)


histogram1 <- histogram_c_chains(second_c_chains_, 2, k2, K2, nclass = 20)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(gamma[1])) + ylab("density")

load("propensity.doublemcmc.RData")
hist_mcmc <- hist(double_mcmc[,2], breaks = histogram1$breaks, plot = F)
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red", alpha = 0.5)
g1

ggsave(filename = "propensity.histogram.gamma1.pdf", plot = g1, width = 5, height = 5)
