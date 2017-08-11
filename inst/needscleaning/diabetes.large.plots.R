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
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

# or load bigger data
load(file = "diabetes.large.more.cchains.RData")
# load(file = "diabetes.cchains_augmented.RData")
X <- scale(diabetes$x2)
p <- ncol(X)
head(df)

df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
df$mean <- (df$sum_est/df$nsamples)


# g <- ggplot(df, aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
# label.df <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(3e-5, p))
# g <- g + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6, position = position_dodge(width = 0.1))
# g <- g + scale_x_log10() + theme(legend.position = "none")
# g
# # g + geom_vline(xintercept = 1e-1)

g <- ggplot(df, aes(x = lambda, y = mean, group = component)) + geom_line() + geom_point()
g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
g <- g + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)) + theme(legend.position = "none")
g <- g + xlab(expression(lambda)) + ylab(TeX("$E_{\\lambda} \\lbrack \\beta | Y,X \\rbrack $"))
# label.df <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(3e-5, p))
components <- c(14, 50, 52, 53, 15, 60, 64)
label.df <- data.frame(component = components,  y = as.numeric(lm(Y ~ X)$coef)[1+components], x = rep(3e-5, length(components)))
g <- g + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6)
g

ggsave(filename = "diabetes.lassopath.p64.pdf", plot = g, width = 10, height = 7)

## graph as a function of L1 norm
# maxL1norm <- sum(abs(lm(Y ~ X)$coef))
# normdf <- df %>% group_by(ilambda, lambda) %>% summarise(relL1norm = sum(abs(mean)) / maxL1norm) %>% ungroup()
# df <- merge(df, normdf, by = c("ilambda", "lambda"))
# g <- ggplot(df %>% filter(component <= 10), aes(x = relL1norm, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = relL1norm))
# g <- g + xlab(TeX("$||\\beta_{\\lambda}||_1 /\\max_{\\lambda}(||\\beta_{\\lambda}||_1)$")) + ylab(TeX("$E(\\beta_{\\lambda}|X,Y)$"))
# g <- g + theme(legend.position = "none")
# label.df <- data.frame(component = 1:10,  y = as.numeric(lm(Y ~ X)$coef)[2:(10+1)], x = rep(1.05, 10))
# # label.df <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(1.05, p))
# # label.df$y[3] <- label.df$y[3] + 0.03
# # label.df$y[10] <- label.df$y[10] - 0.03
# # label.df$y[1] <- label.df$y[1] - 0.035
# g <- g + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6, position = position_dodge(width = 0.1))
# g <- g + scale_color_tableau() #+ scale_x_log10()
# g
# ggsave(filename = "blasso.paths.pdf", plot = g, width = 8, height = 7)


df.p64 <- df

gmeeting <- ggplot(df.p64, aes(x = lambda, y = k-1)) + geom_line() + geom_point() + scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4))
gmeeting <- gmeeting + xlab(TeX("$\\lambda$")) + ylab("90% quantile of meeting times")
gmeeting <- gmeeting + scale_y_log10(limits = c(1, 1e4), breaks = c(5, 10, 30, 100, 1000, 10000))
gmeeting
ggsave(filename = "diabetes.meetingtimes.p64.pdf", plot = gmeeting, width = 10, height = 7)

load("diabetes.cchains.RData")
gmeeting <- gmeeting + geom_line(data=df, aes(x = lambda, y = k-1)) + geom_point(data=df, aes(x = lambda, y = k-1))
gmeeting
label.df <- data.frame(x = c(1e2, 1e3),  y = c(5, 1000), lab = c("p = 10", "p = 64"))
gmeeting <- gmeeting + geom_label(data = label.df, aes(x = x, y = y, label = lab), size = 8)
gmeeting
ggsave(filename = "diabetes.meetingtimes.both.pdf", plot = gmeeting, width = 10, height = 7)

# gmeeting <- gmeeting + scale_y_log10(limits = c(1, 100), breaks = c(1, 3, 10, 30, 100))
# gmeeting
# ggsave(filename = "diabetes.meetingtimes.p10.pdf", plot = gmeeting, width = 10, height = 7)

# df %>% tail
# gmeeting <- ggplot(df, aes(x = lambda, y = k)) + geom_line() + geom_point() + scale_x_log10(breaks = c(1e-4, 1e-2, 1e0, 1e2, 1e4))
# gmeeting <- gmeeting + xlab(TeX("$\\lambda$")) + ylab("k chosen by heuristics")
# gmeeting <- gmeeting +  scale_y_log10(breaks = c(10,100,1000,10000))
# gmeeting + geom_line(aes(y = meantau)) + geom_line(aes(y = mediantau))
# ggsave(filename = "blasso.meeting.pdf", plot = gmeeting, width = 8, height = 7)
# ggplot(df, aes(x = lambda, y = sd/abs(mean), group = component, colour = factor(component))) + geom_line() + geom_point() + scale_x_log10()

# lars_ <- lars::lars(X, Y, trace=TRUE)
# plot(lars_)


### with more nsamples
load(file = "diabetes.large.more.more.cchains.RData")
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



# g <- ggplot(df, aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
# label.df <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(3e-5, p))
# g <- g + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6, position = position_dodge(width = 0.1))
# g <- g + scale_x_log10() + theme(legend.position = "none") + xlab(expression(lambda)) + ylab(TeX("$E(\\beta_{\\lambda}|X,Y)$"))
# g
#
#
# maxL1norm <- sum(abs(lm(Y ~ X)$coef))
# normdf <- df %>% group_by(ilambda, lambda) %>% summarise(relL1norm = sum(abs(mean)) / maxL1norm) %>% ungroup()
# df <- merge(df, normdf, by = c("ilambda", "lambda"))
# g <- ggplot(df %>% filter(component <= 10), aes(x = relL1norm, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = relL1norm))
# g <- g + xlab(TeX("$||\\beta_{\\lambda}||_1 /\\max_{\\lambda}(||\\beta_{\\lambda}||_1)$")) + ylab(TeX("$E(\\beta_{\\lambda}|X,Y)$"))
# g <- g + theme(legend.position = "none")
# label.df <- data.frame(component = 1:10,  y = as.numeric(lm(Y ~ X)$coef)[2:(10+1)], x = rep(1.05, 10))
# # label.df <- data.frame(component = 1:p,  y = as.numeric(lm(Y ~ X)$coef)[2:(p+1)], x = rep(1.05, p))
# # label.df$y[3] <- label.df$y[3] + 0.03
# # label.df$y[10] <- label.df$y[10] - 0.03
# # label.df$y[1] <- label.df$y[1] - 0.035
# g <- g + geom_label(data = label.df, aes(x = x, y = y, label = component), size = 6, position = position_dodge(width = 0.1))
# g <- g + scale_color_tableau() #+ scale_x_log10()
# g

