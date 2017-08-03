# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
setwd("~/Dropbox/PolyaGammaResults/mixture/")

## target distribution
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

load(file = "mixture.c_chains.RData")
nsamples <- length(c_chains_)
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime, nclass = 100)
x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g
# ggsave(filename = "mixture.meetingtime.pdf", plot = g, width = 7, height = 7)



##
k <- 50
K <- 200
histogram <- histogram_c_chains(c_chains_continued_, 1, k, K, nclass = 100)
g <- plot_histogram(histogram, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "mixture.histogram1.pdf", plot = g, width = 7, height = 7)


meetingtime2 <- sapply(c_chains_2, function(x) x$meetingtime)
sum(meetingtime2>100)
# if we are unlucky, none of the large variance is seeable from first 1000
# e.g.
indices <- sample(which(meetingtime2<100), size = 1000)
summary(meetingtime2[indices])

histogram2 <- histogram_c_chains(c_chains_continued_2[indices], 1, k, K, nclass = 100)
g <- plot_histogram(histogram2, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "mixture.histogram2.pdf", plot = g, width = 7, height = 7)

summary(meetingtime2)
histogram2 <- histogram_c_chains(c_chains_continued_2, 1, k, K, nclass = 100)
g <- plot_histogram(histogram2, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "mixture.histogram3.pdf", plot = g, width = 7, height = 7)

