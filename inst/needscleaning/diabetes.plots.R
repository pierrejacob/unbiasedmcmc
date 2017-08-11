# load packages
library(debiasedmcmc)
library(ggthemes)
library(lars)
library(latex2exp)
setmytheme()
rm(list = ls())
set.seed(21)
# registerDoParallel(cores = detectCores())
registerDoParallel(cores = 1)
setwd("~/Dropbox/PolyaGammaResults/reproduce/")

data(diabetes)
X <- scale(diabetes$x)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)
load(file = "diabetes.cchains.RData")

# or load bigger data
# load(file = "diabetes.large.cchains.alt.RData")
# load(file = "diabetes.cchains_augmented.RData")
# X <- scale(diabetes$x2)
# p <- ncol(X)
# head(df)

# df$signal <- df$component <= 5
df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
df$mean <- (df$sum_est/df$nsamples)
g <- ggplot(df, aes(x = lambda, y = mean, group = component)) + geom_line() + geom_point()
g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
g <- g + scale_x_log10(breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3)) + theme(legend.position = "none")
g <- g + xlab(expression(lambda)) + ylab(TeX("$E_{\\lambda} \\lbrack \\beta | Y,X \\rbrack $"))
label.df <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(5e-3, p))
label.df$y[3] <- label.df$y[3] + 0.04
label.df$y[10] <- label.df$y[10] - 0.04
label.df$y[1] <- label.df$y[1] - 0.045
g <- g + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6)
g
ggsave(filename = "diabetes.lassopath.p10.pdf", plot = g, width = 10, height = 7)

## graph as a function of L1 norm
maxL1norm <- sum(abs(lm(Y ~ X)$coef))
normdf <- df %>% group_by(ilambda, lambda) %>% summarise(relL1norm = sum(abs(mean)) / maxL1norm) %>% ungroup()
df <- merge(df, normdf, by = c("ilambda", "lambda"))
gnorm <- ggplot(df, aes(x = relL1norm, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
gnorm <- gnorm + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = relL1norm))
gnorm <- gnorm + xlab(TeX("$||\\beta_{\\lambda}||_1 /\\max_{\\lambda}(||\\beta_{\\lambda}||_1)$")) + ylab(TeX("$E(\\beta_{\\lambda}|X,Y)$"))
gnorm <- gnorm + theme(legend.position = "none")
label.df <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(1.05, p))
label.df$y[3] <- label.df$y[3] + 0.03
label.df$y[10] <- label.df$y[10] - 0.03
label.df$y[1] <- label.df$y[1] - 0.035
gnorm <- gnorm + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6)
gnorm <- gnorm + scale_color_tableau() # + scale_x_log10()
gnorm

# ggsave(filename = "blasso.paths.pdf", plot = g, width = 8, height = 7)
gmeeting <- ggplot(df, aes(x = lambda, y = k-1)) + geom_line() + geom_point() + scale_x_log10(breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3))
gmeeting <- gmeeting + xlab(TeX("$\\lambda$")) + ylab("90% quantile of meeting times")
gmeeting <- gmeeting + scale_y_log10(limits = c(1, 100), breaks = c(1, 3, 10, 30, 100))
gmeeting
ggsave(filename = "diabetes.meetingtimes.p10.pdf", plot = gmeeting, width = 10, height = 7)


# gmeeting <- gmeeting + geom_line(aes(y = mediantau), linetype = 3) + geom_line(aes(y = meantau), linetype = 2)


# ggsave(filename = "blasso.meeting.pdf", plot = gmeeting, width = 8, height = 7)
# ggplot(df, aes(x = lambda, y = sd/abs(mean), group = component, colour = factor(component))) + geom_line() + geom_point() + scale_x_log10()
# lars_ <- lars::lars(X, Y, trace=TRUE)
# plot(lars_)
