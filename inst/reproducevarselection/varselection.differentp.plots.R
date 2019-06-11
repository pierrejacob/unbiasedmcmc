library(unbiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
#
library(doRNG)
library(doParallel)
library(tidyr)
library(viridis)
library(dplyr)
registerDoParallel(cores = detectCores()-2)
#

# later put n = 500 and p = 1000
n <- 500
SNR <- 1
#
# load data
load(paste0("varselection.dataSNR", SNR, ".RData"))
X_full <- X
Y_full <- Y
kappa <- 2
# subset data to desired size
#
s0 <- 100
proportion_singleflip <- 0.5

nrep <- 1000
nps <- 5
ps <- c(100, 250, 500, 750, 1000)
df.ue <- data.frame()
meetingsfilepath <- paste0("varselection.meetings.differentp.n", n, ".RData")

load(meetingsfilepath)
df.summary <- df.ue %>% group_by(ip,p) %>% summarise(m = mean(meetings)) %>% mutate(scaledm = m/p)

g <- ggplot(df.ue, aes(x = p, y = meetings/p, group = p)) + geom_violin() + ylab("meeting times / p")
g <- g + scale_x_continuous(breaks = ps)
g

ggsave(filename = "varselection.meetings.differentp.pdf", plot = g, width = 8, height = 6)
