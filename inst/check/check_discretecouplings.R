# This script plays with coupling of distributions defined on discrete spaces
# load packages
library(debiasedmcmc)
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
#
p <- 40
s <- 20
selection <- rep(0, p)
selection[sample(1:p, s, replace = F)] <- (runif(s) < 0.5)
selection
debiasedmcmc:::sample_pair01(selection)
# this is meant to sample one zero and one one uniformly
nrep <- 10000

test <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  debiasedmcmc:::sample_pair01(selection)
}
#
table(test[,1])  / nrep
which(selection == 0)
1/sum(selection==0)
#

selection1 <- rep(0, p)
selection1[sample(1:p, s, replace = F)] <- (runif(s) < 0.5)
selection2 <- rep(0, p)
selection2[sample(1:p, s, replace = F)] <- (runif(s) < 0.5)


test <- foreach(irep = 1:10000, .combine = rbind) %dorng% {
  coupled_pairs01(selection1, selection2)
}

#
table(test[,1]) / nrow(test)
1/(p - sum(selection1))
table(test[,2]) / nrow(test)
1/(p - sum(selection2))
table(test[,3]) / nrow(test)
1/sum(selection1)
table(test[,4]) / nrow(test)
1/sum(selection2)
#

