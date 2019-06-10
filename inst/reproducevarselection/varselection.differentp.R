library(debiasedmcmc)
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

for (ip in seq_along(ps)){
  print(ip)
  p <- ps[ip]
  Y <- Y_full[1:n]
  X <- X_full[1:n,1:p]
  Y2 <- (t(Y) %*% Y)[1,1]
  g <- p^3
  # load model
  vs <- get_variableselection(Y,X,g,kappa,s0,proportion_singleflip)
  prior <- vs$prior
  marginal_likelihood <- vs$marginal_likelihood
  rinit <- vs$rinit
  single_kernel <- vs$single_kernel
  coupled_kernel <- vs$coupled_kernel
  unbiasedestimator <- vs$unbiasedestimator
  # Unbiased MCMC
  ues <- foreach(i = 1:nrep) %dorng% {
    unbiasedestimator(single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 0)
  }
  meetings <- sapply(ues, function(x) x$meetingtime)
  df.ue <- rbind(df.ue, data.frame(irep = 1:nrep, ip = rep(ip, nrep), p = rep(p, nrep), meetings = meetings))
  save(df.ue, file = meetingsfilepath)
}

load(meetingsfilepath)
df.summary <- df.ue %>% group_by(ip,p) %>% summarise(m = mean(meetings)) %>% mutate(scaledm = m/p)

g <- ggplot(df.ue, aes(x = p, y = meetings/p, group = p)) + geom_violin() + ylab("meeting times / p")
g <- g + scale_x_continuous(breaks = ps)
g

ggsave(filename = "varselection.meetings.differentp.pdf", plot = g, width = 8, height = 6)


