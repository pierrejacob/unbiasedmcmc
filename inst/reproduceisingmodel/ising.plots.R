### This script creates plots
library(debiasedmcmc)
library(viridis)
library(dplyr)
rm(list = ls())
set.seed(21)
setmytheme()
#

## Load meeting times associated with single-site Gibbs sampler
load(file = "ising.singlesite.meetings.RData")
tail(singlesite.meetings.df)
unique(singlesite.meetings.df$beta)
g <- ggplot(singlesite.meetings.df %>% group_by(beta) %>% summarise(m = mean(meeting)),
            aes(x = beta, y = m)) + geom_line() + geom_point() + scale_y_log10(limits = c(1e1, 1e7), breaks = c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7))
g <- g + xlab(expression(theta)) + ylab("average meeting time")
g

ggsave(filename = "ising.singlesite.meetings.pdf", plot = g, width = 8, height = 6)

## Load meetings associated with parallel tempering

nchains_values <- c(4,8,12,16,24,32)
nchains.df <- data.frame()
for (nchains in nchains_values){
  load(paste0("ising.swapN", nchains, ".meetings.RData"))
  nchains.df <- rbind(nchains.df, data.frame(nchains = nchains, mean = mean(meetings),
                                             q90 = as.numeric(quantile(meetings, probs = 0.9)),
                                             max = max(meetings)))
}
nchains.df

#
g <- ggplot(nchains.df, aes(x = nchains, y = mean))
g <- g + geom_line() + geom_point() +  scale_y_log10(limits = c(1e1, 1e7), breaks = c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7))
g <- g + xlab("# of chains") + ylab("average meeting time") + scale_x_continuous(breaks = nchains_values)
g
ggsave(filename = "ising.parallel.meetings.pdf", plot = g, width = 8, height = 6)

### load results from cluster
# datafiles <- list.files(path = "output/", pattern = "5000")
datafiles <- list.files(path = "output/", pattern = "1e+")
nchains <- 16
betas <- seq(from = 0.3, to = 0.55, length.out = nchains)
datafiles[1]
nchains <- 16
df_full <- data.frame()
for (ifile in seq_along(datafiles)){
  load(file = paste0("output/", datafiles[ifile]))
  uestimators <- t(sapply(results_, function(x) x$uestimator))
  starttime <- c(0, cumsum(durations_))[1:(length(durations_))]
  endtime <- cumsum(durations_)
  df_full <- rbind(df_full, data.frame(rep = rep(ifile, nsamples_),
                             isample = 1:nsamples_,
                             durations = durations_,
                             starttime = starttime,
                             endtime = endtime,
                             meetings = sapply(results_, function(x) x$meetingtime),
                             uestimators = uestimators))
}

max(df_full$endtime)/3600
# subset df by completion time in seconds
df <- df_full %>% filter(endtime < 30*60)


g <- ggplot(df, aes(y = rep, yend = rep, x = starttime/60+0.01, xend = endtime/60-0.01, colour = isample)) + geom_segment(lineend = "round", alpha = 1)
g <- g + geom_point() +  scale_color_viridis(name = "sample index:", discrete = FALSE)
g <- g + theme(legend.position = "none") + xlab("time (minutes)") + ylab("processor")  + geom_vline(xintercept = max(df$endtime)/60, linetype = 2)
g
ggsave(filename = "ising.chronology.png", plot = g, width = 8, height = 6)

# sort(unique((df %>% group_by(rep) %>% summarise(max = max(isample)))$max))
nrow(df)
# mean(df$durations)/60

# mean((df %>% group_by(rep) %>% summarise(max = max(isample)))$max)
ggplot(df, aes(x = meetings, y = durations/60, colour = rep)) + geom_point() + ylab("minutes") +
  scale_color_viridis(discrete = FALSE) + theme(legend.position = "none")

# to get meeting times for a histogram without bias,
# we can simply get the first X samples of each processors

meetings <- (df_full %>% filter(isample <= 10) %>% select(meetings))$meetings
qplot(x = meetings, geom = "blank") + geom_histogram(aes(y = ..density..))

cat("average meeting time:", mean(meetings), "\n")
# mean(meetings < 5e4)
# quantile(df$meetings, probs = c(0.9, 0.95, 0.99))
# compute average per processor, for each beta
df.summary <- df %>% group_by(rep) %>% summarise_at(.vars = names(.)[7:22],
                                      .funs = c(mean = "mean"))

# compute average across processors
estimators <- as.numeric(colMeans(df.summary[,2:17]))
sd_estimators <- as.numeric(apply(df.summary[,2:17], 2, sd))
# the above is \hat{\sigma}_1(P,t) in the article

df.plot <- data.frame(betas = betas, estimators = estimators, sd_estimators = sd_estimators)
df.plot$shat.eq34 <- sd_estimators/sqrt(nrow(df.summary))
# confidence intervals valid when P goes to infinity, fixed t
# i.e. Eq. (3.4)
ggplot(df.plot, aes(x = betas, y = estimators)) + geom_line() +
  geom_errorbar(aes(ymin = estimators - 2*shat.eq34, ymax = estimators + 2*shat.eq34))

# but we can compute the standard deviation differently
# for it to be valid for fixed P, t goes to infinity, Eq. (3.5)

mean_durations <- mean((df %>% group_by(rep) %>% summarise(m = mean(durations)))$m)
averagesquare <- as.numeric(colMeans(df %>% group_by(rep) %>% summarise_at(.vars = names(.)[7:22], .funs = function(x) mean(x^2)))[2:17])
sigma2hat <- mean_durations * (averagesquare - estimators^2)
# so compute width of estimators as
max(df$endtime)
df.plot$shat.eq35 <- sqrt(sigma2hat) / sqrt(500 * max(df$endtime))

g <- ggplot(df.plot, aes(x = betas, y = estimators)) + geom_line() + xlab(expression(theta)) + geom_point()
  # geom_errorbar(aes(ymin = estimators - 2*shat.eq35, ymax = estimators + 2*shat.eq35)) +
g <- g + scale_x_continuous(breaks = seq(from = 0.3, to = 0.55, by = 0.05)) + ylab("natural stat.")
g

# comparison of standard deviation, Eqs. (3.4) and (3.5)
g2 <- ggplot(df.plot, aes(x = betas, y = shat.eq34)) + geom_line() + geom_point() # + geom_line(aes(y = shat.eq35), colour = "red")
# g <- g + geom_vline(xintercept = 1/2.29, linetype = 2)
g2 <- g2 + xlab(expression(theta)) + scale_x_continuous(breaks = seq(from = 0.3, to = 0.55, by = 0.05))
g2 <- g2 + ylab("standard error")
g2

g12 <- gridExtra::grid.arrange(g, g2, ncol = 1)
ggsave(filename = "ising.estimates.pdf", plot = g12, width = 8, height = 6)
# the max seems to correspond roughly to critical temperature?



