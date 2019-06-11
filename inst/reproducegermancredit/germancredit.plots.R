
# load packages
library(unbiasedmcmc)
library(dplyr)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())
data(germancredit)
load("germancredit.tuning.RData")

#ndata <- 10


# Plot histogram of meetingtimes
nsamples <- meetingtime.list$nsamples
meetingtimes <- sapply(X = meetingtime.list$meetings_1, FUN = function(x) x$meetingtime)
#upper quantiles: 90% = 81, 95% = 96, 99% = 135. In the paper we consider 110.

mt_df = data.frame(meetingtime=meetingtimes)
head(mt_df)
#mt_df <- data.frame(meetingtime=sapply(c_chains, function(x) x$meetingtime))
g <- ggplot(data=mt_df,aes(x=meetingtime))+geom_histogram(aes(y = ..density..), binwidth=10)
g <- g + xlab("meeting time")
g
ggsave(filename="pgg.meetingtime.pdf", plot=g, width=7, height=7)



# Plot distances and number met for a particular run
g3 <- ggplot(df_w,aes(x=iter,y=nmet_w))+geom_line() + xlab('iteration') + ylab("count")
g3
ggsave(filename = "pgg.nmetw.pdf", plot = g3, width = 7, height = 7)

g2 <- ggplot(df_w,aes(x=iter,y=dist_w))+geom_line() + xlab('iteration') + ylab('distance') #+ ylab(TeX("$||w_1-w_2||_2$"))
g2
ggsave(filename = "pgg.distw.pdf", plot = g2, width = 7, height = 7)

g1 <- ggplot(df_beta,aes(x=iter,y=dist_beta))+geom_line() + xlab('iteration') + ylab('distance') #ylab(TeX("$||\\beta_1-\\beta_2||_2$"))
g1
ggsave(filename = "pgg.distbeta.pdf", plot = g1, width = 7, height = 7)




# Plot efficiency as a function of k
tuning.k.list$nsamples
g <- qplot(x = tuning.k.list$ks, y = 1/(tuning.k.list$cost * tuning.k.list$v), geom = "line") #+ scale_y_continuous(breaks = c(0.2, 0.3, 0.4))
g <- g + geom_point() + ylab('asymptotic efficiency') + xlab("k")
g
ggsave(filename = "pgg.asympt_eff.pdf", plot = g, width = 7, height = 7)






# Plot histogram

load(file = "germancredit.c_chain.RData")
load(file = "germancredit.mcmc.RData")

idx1 <- which('Instalment.per.cent' == colnames(X))
idx2 <- which('Duration.in.Current.address' == colnames(X))
idx <- idx1

nclass <- 27 #27 for idx1, 24 for idx2
rng <- range(find_breaks(c_chains_, idx, nclass, k, m, lag = 1))
breaks = seq(rng[1],rng[2],length=nclass)

histogram1 <- histogram_c_chains(c_chains_, idx, k, m, nclass = nclass)

niterations <- nrow(chain)
hist_mcmc <- hist(chain[1000:niterations, idx], breaks = histogram1$breaks, plot = FALSE)

g1 <- plot_histogram(histogram1, with_bar = TRUE) + xlab(expression(beta)) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
g1
ggsave(filename = "pgg.histogram1.pdf", plot = g1, width = 7, height = 7)

library(coda)
testfunction <- function(x) x[idx1]
mcmc_ <- apply(X = chain[1000:nrow(chain),,drop=F], MARGIN = 1, FUN = testfunction)
mcmcvar <- spectrum0(mcmc_)$spec
mcmcvar

nsamples <- length(c_chains_)
k <- 110
m <- 1100
taus <- sapply(c_chains_, FUN = function(x) x$meetingtime)
estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_[[irep]], h = testfunction, k = k, m = m)
}
v <- var(unlist(estimators))
c <- mean(2 * taus + pmax(1, m + 1 - taus))
inef <- c * v
# loss of efficiency for this choice of k and m
inef / mcmcvar
#

