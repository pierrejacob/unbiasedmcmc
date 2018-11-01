# load packages
library(debiasedmcmc)
library(latex2exp)
library(dplyr)
setmytheme()
rm(list = ls())
set.seed(21)

data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

# lambdas <- 10^(seq(from = -2, to = 3, length.out = 25))

## meeting times as a function of lambda
load("bayesianlasso.meetings.RData")
head(df)
meetingquantile <- df %>% group_by(ilambda) %>% summarise(lambda = mean(lambda), mean = mean(meetingtime), q99 = quantile(meetingtime, probs = 0.99))
meetingquantile %>% head
g <- ggplot(meetingquantile, aes(x = lambda, y = mean)) + geom_point()
g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4))
g <- g + scale_y_log10(breaks = c(10, 100, 1000, 5500))
g <- g + xlab(expression(lambda)) + ylab("average meeting time")
g
ggsave(filename = "bayesianlasso.meetings.pdf", plot = g, width = 8, height = 6)

## effective sample sizes as a function of lambda
load("bayesianlasso.mcmc.RData")
g <- ggplot(df, aes(x = lambda, y = ess / (nmcmc - burnin + 1))) + geom_point()
g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4))
g <- g + xlab(expression(lambda)) + ylab("effective sample sizes")
g
ggsave(filename = "bayesianlasso.ess.pdf", plot = g, width = 8, height = 6)

## posterior mean of beta as function of lambda
load("bayesianlasso.estimators.RData")
summary.df <- df %>% group_by(ilambda, lambda, component) %>% summarise(m = mean(estimator), sd = sd(estimator), nsamples = n())
summary.df$component %>% unique
summary.df$nsamples %>% unique

load("bayesianlasso.extraestimators.RData")
summary.extra.df <- rbind(df, extra.df) %>% group_by(ilambda, lambda, component) %>% summarise(m = mean(estimator), sd = sd(estimator), nsamples = n())

# g <- ggplot(summary.df %>% filter(component <= 64), aes(x = lambda, y = sd/sqrt(nsamples), group = component)) + geom_line() + geom_point()
# g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)) + theme(legend.position = "none")
# g <- g + ylim(0,0.4)
# g
#
# g <- ggplot(summary.extra.df %>% filter(component <= 64), aes(x = lambda, y = sd/sqrt(nsamples), group = component)) + geom_line() + geom_point()
# g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)) + theme(legend.position = "none")
# g <- g + ylim(0,0.4)
# g

g <- ggplot(summary.df %>% filter(component <= 64), aes(x = lambda, y = m, group = component)) + geom_line() + geom_point()
g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)) + theme(legend.position = "none")
g <- g + geom_segment(aes(y = m - 2*sd / sqrt(nsamples), yend = m + 2*sd / sqrt(nsamples), xend = lambda))
g <- g + xlab(expression(lambda)) + ylab(TeX("$E_{\\lambda} \\lbrack \\beta | Y,X \\rbrack $"))
g <- g + ylim(-5,4)
g
ggsave(filename = "bayesianlasso.paths.pdf", plot = g, width = 8, height = 6)

g <- ggplot(summary.extra.df %>% filter(component <= 64), aes(x = lambda, y = m, group = component)) + geom_line() + geom_point()
g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)) + theme(legend.position = "none")
g <- g + geom_segment(aes(y = m - 2*sd / sqrt(nsamples), yend = m + 2*sd / sqrt(nsamples), xend = lambda))
g <- g + xlab(expression(lambda)) + ylab(TeX("$E_{\\lambda} \\lbrack \\beta | Y,X \\rbrack $"))
g <- g + ylim(-5,4)
g
ggsave(filename = "bayesianlasso.refinedpaths.pdf", plot = g, width = 8, height = 6)

