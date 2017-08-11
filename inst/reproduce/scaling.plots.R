# load packages
library(debiasedmcmc)
library(ggthemes)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#

load(file = "dimension.scaling.gibbs.itercst.RData")
dimensions <- c(1, 50, 100, 200, 300)
iterates <- c(1, 2, 5)
nsamples <- 1000
iterates_labels <- c("1 step", "2 steps", "5 steps")
df$iterate <- factor(df$iterate, levels = iterates, labels = iterates_labels)
g <- ggplot(df, aes(x = d, y = mean_time, group = iterate, colour = factor(iterate))) + geom_line() + geom_point() + ylab("average meeting time")
# g <- ggplot(df, aes(x = d, y = mean_time)) + geom_line() + scale_x_continuous(breaks = c(1,25, 50,75, 100)) + ylab("average meeting time")
# g <- g + scale_y_continuous(breaks = c(1, 10, 20), limits = c(0,50))
g <- g + scale_color_colorblind() + scale_x_continuous(breaks = dimensions)
g <- g + geom_label(data = data.frame(x = c(250,250,250), y = c(57, 28, 12), iterate = iterates_labels),
                    aes(x = x, y = y, colour = iterate, label = iterate), size = 8) + theme(legend.position = "none")
g <- g + xlab("dimension")
g
# g + scale_x_log10() + scale_y_log10()

ggsave(filename = "dimension.scaling.gibbs.pdf", plot = g, width = 8, height = 7)

load(file = "dimension.scaling.rwmh.RData")
df.scaling1 <- df
# g <- ggplot(df, aes(x = d, y = mean_time)) + geom_line() + scale_x_continuous(breaks = seq(from = 1, to = 15, by = 2)) + ylab("average meeting time")
# g <- g + scale_y_log10() + geom_point()
# g
# ggsave(filename = "dimension.scaling.rwmh1.pdf", plot = g, width = 5, height = 5)

load(file = "dimension.scaling.rwmh.optimal.RData")
df <- df %>% filter(is.finite(max_time))
g <- ggplot(df, aes(x = d, y = mean_time, colour = "scaling 1")) + geom_line() + scale_x_continuous(breaks = seq(from = 1, to = 15, by = 2)) + ylab("average meeting time")
g <- g + scale_y_log10(breaks = c(10, 100, 1000, 10000)) + geom_point() + xlab("dimension")
g <- g + geom_line(data = df.scaling1, aes(colour = "scaling 2")) + geom_point(data = df.scaling1, aes(colour = "scaling 2"))
g <- g + scale_color_colorblind()
g <- g + geom_label(data = data.frame(x = c(5,12), y = c(7000, 50), scaling = c("scaling 1", "scaling 2")),
                    aes(x = x, y = y, colour = scaling, label = scaling), size = 8) + theme(legend.position = "none")
g

ggsave(filename = "dimension.scaling.rwmh.pdf", plot = g, width = 8, height = 7)

