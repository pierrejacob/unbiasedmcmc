# load packages
library(debiasedmcmc)
library(dplyr)
library(tidyr)
library(ggthemes)
setmytheme()
rm(list = ls())
set.seed(21)
#

load(file = "scalingdimension.rwmh.maxcoupling.RData")
df
g <- ggplot(df %>% filter(target_type == "dense"), aes(x = d, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
g <- g + scale_y_log10(breaks = c(10,100,1e3, 1e4)) + xlab("dimension")
g
ggsave(filename = "scalingdimension.rwmh.maxcoupling.pdf", plot = g, width = 8, height = 6)

load(file = "scalingdimension.rwmh.reflmaxcoupling.RData")
g <- ggplot(df, aes(x = d, y = mean_time, group = init_type, linetype = init_type)) + geom_line() + ylab("average meeting time")
g <- g + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
g <- g + xlab("dimension") + scale_linetype("initialization:")
g
ggsave(filename = "scalingdimension.rwmh.reflmaxcoupling.pdf", plot = g, width = 8, height = 6)


# g <- ggplot(df %>% filter(target_type == "dense", init_type == "offset"), aes(x = d, y = mean_time, group = init_type)) + geom_line() + ylab("average meeting time")
# g <- g + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
# g <- g + xlab("dimension")
# g
# ggsave(filename = "scalingdimension.rwmh.reflmaxcoupling.offsetinit.pdf", plot = g, width = 5, height = 5)
#
# g <- ggplot(df %>% filter(target_type == "dense", init_type == "target"), aes(x = d, y = mean_time, group = init_type)) + geom_line() + ylab("average meeting time")
# g <- g + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
# g <- g + xlab("dimension")
# g
# ggsave(filename = "scalingdimension.rwmh.reflmaxcoupling.targetinit.pdf", plot = g, width = 5, height = 5)

# g <- ggplot(df %>% filter(target_type == "dense"), aes(x = d, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g +  scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
# g <- g + xlab("dimension")
# g
# ggsave(filename = "scalingdimension.rwmh.reflmaxcoupling.pdf", plot = g, width = 8, height = 6)


load(file = "scalingdimension.gibbs.sparse.RData")
df.sparse <- df
load("scalingdimension.gibbs.dense.RData")
df.dense <- df

g <- ggplot(df.sparse, aes(x = d, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.sparse$d))) + xlab("dimension")
g
ggsave(filename = "scalingdimension.gibbs.sparse.pdf", plot = g, width = 8, height = 6)

g <- ggplot(df.dense, aes(x = d, y = median_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("median meeting time")
g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.dense$d))) + xlab("dimension")
# g <- g + scale_y_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5))
g
ggsave(filename = "scalingdimension.gibbs.dense.pdf", plot = g, width = 8, height = 6)

load("scalingdimension.wishart.hmc.meetings.RData")
df.summary <- df %>% group_by(dimension, init_type) %>% summarise(mean_time = mean(meetingtimes))
g <- ggplot(df.summary, aes(x = dimension, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.summary$dimension))) + xlab("dimension") + ylim(0, 60)
g

# ggsave(filename = "scalingdimension.hmc.dense.pdf", plot = g, width = 8, height = 6)
