# load packages
library(debiasedmcmc)
library(ggthemes)
library(latex2exp)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())

data(diabetes)
X <- scale(diabetes$x)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)
load(file = "diabetes.cchains.p10.RData")
df.p10 <- df

df.p10$sd <- sqrt(df.p10$sum_squareest/df.p10$nsamples - (df.p10$sum_est/df.p10$nsamples)^2)
df.p10$mean <- (df.p10$sum_est/df.p10$nsamples)
g <- ggplot(df.p10, aes(x = lambda, y = mean, group = component)) + geom_line() + geom_point()
g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
g <- g + scale_x_log10(breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3)) + theme(legend.position = "none")
g <- g + xlab(expression(lambda)) + ylab(TeX("$E_{\\lambda} \\lbrack \\beta | Y,X \\rbrack $"))
label.df.p10 <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(5e-3, p))
label.df.p10$y[3] <- label.df.p10$y[3] + 0.04
label.df.p10$y[10] <- label.df.p10$y[10] - 0.04
label.df.p10$y[1] <- label.df.p10$y[1] - 0.045
g <- g + geom_label(data = label.df.p10, aes(x = x, y = y, label = component), size = 6)
g
ggsave(filename = "diabetes.lassopath.p10.pdf", plot = g, width = 10, height = 7)

gmeeting <- ggplot(df.p10, aes(x = lambda, y = k-1)) + geom_line() + geom_point() + scale_x_log10(breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3))
gmeeting <- gmeeting + xlab(TeX("$\\lambda$")) + ylab("90% quantile of meeting times")
gmeeting <- gmeeting + scale_y_log10(limits = c(1, 100), breaks = c(1, 3, 10, 30, 100))
gmeeting

load(file = "diabetes.cchains.p64.RData")
df.p64 <- df
data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

df.p64$sd <- sqrt(df.p64$sum_squareest/df.p64$nsamples - (df.p64$sum_est/df.p64$nsamples)^2)
df.p64$mean <- (df.p64$sum_est/df.p64$nsamples)

g <- ggplot(df.p64, aes(x = lambda, y = mean, group = component)) + geom_line() + geom_point()
g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)) + theme(legend.position = "none")
g <- g + xlab(expression(lambda)) + ylab(TeX("$E_{\\lambda} \\lbrack \\beta | Y,X \\rbrack $"))
components <- c(14, 50, 52, 53, 15, 60, 64)
label.df.p64 <- data.frame(component = components,  y = as.numeric(lm(Y ~ X)$coef)[1+components], x = rep(3e-5, length(components)))
g <- g + geom_label(data = label.df.p64, aes(x = x, y = y, label = component), size = 6)
g
ggsave(filename = "diabetes.lassopath.p64.pdf", plot = g, width = 10, height = 7)

gmeeting <- ggplot(df.p64, aes(x = lambda, y = k-1)) + geom_line() + geom_point() + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4))
gmeeting <- gmeeting + xlab(TeX("$\\lambda$")) + ylab("90% quantile of meeting times")
gmeeting <- gmeeting + scale_y_log10(limits = c(1, 1e4), breaks = c(5, 10, 30, 100, 1000, 10000))
gmeeting <- gmeeting + geom_line(data=df.p10, aes(x = lambda, y = k-1)) + geom_point(data=df, aes(x = lambda, y = k-1))
label.df <- data.frame(x = c(1e2, 1e3),  y = c(5, 1000), lab = c("p = 10", "p = 64"))
gmeeting <- gmeeting + geom_label(data = label.df, aes(x = x, y = y, label = lab), size = 8)
gmeeting
ggsave(filename = "diabetes.meetingtimes.both.pdf", plot = gmeeting, width = 10, height = 7)

## Then plot refined estimates

load(file = "diabetes.cchains.p64.refine.RData")
df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
df$mean <- (df$sum_est/df$nsamples)

g <- ggplot(df, aes(x = lambda, y = mean, group = component)) + geom_line() + geom_point()
g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)) + theme(legend.position = "none")
g <- g + xlab(expression(lambda)) + ylab(TeX("$E_{\\lambda} \\lbrack \\beta | Y,X \\rbrack $"))
components <- c(14, 50, 52, 53, 15, 60, 64)
label.df <- data.frame(component = components,  y = as.numeric(lm(Y ~ X)$coef)[1+components], x = rep(3e-5, length(components)))
g <- g + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6)
g

ggsave(filename = "diabetes.lassopath.p64.refined.pdf", plot = g, width = 10, height = 7)
