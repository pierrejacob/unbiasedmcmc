# load packages
library(debiasedmcmc)
library(ggthemes)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#


load(file = "scalingdimension.gibbs.RData")
dimensions <- c(1, 50, 100, 200, 300)
iterates <- c(1, 2, 5)
iterates_labels <- c("1 step", "2 steps", "5 steps")
df$iterate <- factor(df$iterate, levels = iterates, labels = iterates_labels)
g <- ggplot(df, aes(x = d, y = mean_time, group = iterate, colour = factor(iterate))) + geom_line() + geom_point() + ylab("average meeting time")
g <- g + scale_color_colorblind() + scale_x_continuous(breaks = dimensions)
g <- g + geom_label(data = data.frame(x = c(250,250,250), y = c(57, 28, 12), iterate = iterates_labels),
                    aes(x = x, y = y, colour = iterate, label = iterate), size = 8) + theme(legend.position = "none")
g <- g + xlab("dimension")
g
ggsave(filename = "scalingdimension.gibbs.pdf", plot = g, width = 6, height = 5)

load(file = "scalingdimension.rwmh.scaling1.RData")
df.scaling1 <- df
# g <- ggplot(df, aes(x = d, y = mean_time)) + geom_line() + scale_x_continuous(breaks = seq(from = 1, to = 15, by = 2)) + ylab("average meeting time")
# g <- g + scale_y_log10() + geom_point()
# g
# ggsave(filename = "scalingdimension.rwmh1.pdf", plot = g, width = 5, height = 5)

load(file = "scalingdimension.rwmh.scaling2.RData")
df.scaling2 <- df

g <- ggplot(df.scaling1, aes(x = d, y = mean_time, colour = "scaling 1")) + geom_line() + scale_x_continuous(breaks = seq(from = 1, to = 15, by = 2)) + ylab("average meeting time")
g <- g + scale_y_log10(breaks = c(10, 100, 1000, 10000)) + geom_point() + xlab("dimension")
g <- g + geom_line(data = df.scaling2, aes(colour = "scaling 2")) + geom_point(data = df.scaling2, aes(colour = "scaling 2"))
g <- g + scale_color_colorblind()
g <- g + geom_label(data = data.frame(x = c(5,12), y = c(7000, 50), scaling = c("scaling 1", "scaling 2")),
                    aes(x = x, y = y, colour = scaling, label = scaling), size = 8) + theme(legend.position = "none")
g
ggsave(filename = "scalingdimension.rwmh.pdf", plot = g, width = 6, height = 5)

load(file = "scalingdimension.hmc.RData")
df$effectiveTime %>% unique
l <- levels(factor(df$effectiveTime))
df$effectiveTime <- factor(df$effectiveTime, levels = l, labels = c("pi/4", "pi/3", "pi/2"))
# df$effectiveTime <- factor(df$effectiveTime, levels = l, labels = c("pi/8", "pi/6", "pi/3"))
levels(df$effectiveTime)
g <- ggplot(df, aes(x = d, y = mean_time, group = effectiveTime, colour = factor(effectiveTime))) + geom_line() + geom_point() + ylab("average meeting time")
g <- g + scale_color_colorblind() + scale_x_continuous(breaks = dimensions) + theme(legend.position = "none")
g <- g + xlab("dimension") #+ ylim(0, 200)
g <- g + geom_label(data = data.frame(x = c(290, 290, 290), y = c(105, 45, 130), effectiveTime = c("pi/4", "pi/3", "pi/2")),
                    aes(x = x, y = y, colour = effectiveTime, label = effectiveTime), size = 8, parse = T)
# g <- g + geom_label(data = data.frame(x = c(290, 290, 290), y = c(375, 200, 40), effectiveTime = c("pi/8", "pi/6", "pi/3")),
#                     aes(x = x, y = y, colour = effectiveTime, label = effectiveTime), size = 8, parse = T)
g

ggsave(filename = "scalingdimension.hmc.pdf", plot = g, width = 6, height = 5)
