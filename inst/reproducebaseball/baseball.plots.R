# load packages
library(unbiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores())
# this example is taken from Rosenthal "Parallel Computing and Monte Carlo algorithms", 2000.
# The data are baseball player's batting averages, as in Efron and Morris 1975 or Morris 1983
# In particular, the data are listed in Morris 1983, and available as part of Rosenthal's online code:
# http://probability.ca/jeff/comp/james.c

Y <- c(0.395, 0.375, 0.355, 0.334, 0.313, 0.313, 0.291, 0.269, 0.247, 0.247, 0.224, 0.224,
       0.224, 0.224, 0.224, 0.200, 0.175, 0.148)
ndata <- length(Y)

# The data are Y_i for i in {1,...,K}, with here K = 18.
# The model specifies Y_i ~ Normal(theta_i, V), with V known and fixed to 0.00434
V <- 0.00434
# The prior is theta_i ~ Normal(mu, A), with mu ~ flat prior, and A ~ InverseGamma(a,b)
# with density proportional to exp(-b/x) x^{-a-1}; Rosenthal chooses a = -1, b = 2
a <- -1
b <- 2
# To target the posterior distribution, we follow Rosenthal and consider the following Gibbs sampler.
# A given rest: IG(a + (K-1)/2, b + (sum_i=1^K (theta_i - theta_bar)^2)/2)
# mu given rest: Normal(theta_bar, A / K)
# theta_i given rest: Normal( (mu * V + Y_i * A) / (V + A), A * V / (V + A))
# ... where theta_bar is the average of the theta_i, for i in {1,...,K}

load(file = "baseball.tuning.RData")

testfunction <- function(x) x[1]^2

load(file = "baseball.mcmc.RData")
library(coda)
#
mcmc_ <- apply(X = chain[1000:nrow(chain),,drop=F], MARGIN = 1, FUN = testfunction)
mcmcvar <- spectrum0(mcmc_)$spec


nsamples <- length(c_chains_)
taus <- sapply(c_chains_, FUN = function(x) x$meetingtime)
table(taus)
ks <- c(1, 3, 5)
mfactor <- c(1, 5, 10)
ineff.df <- data.frame()
for (ik in 1:length(ks)){
  for (im in 1:length(mfactor)){
    k <- ks[ik]
    m <- mfactor[im] * k
    estimators <-  foreach(irep = 1:nsamples) %dorng% {
      H_bar(c_chains_[[irep]], h = testfunction, k = k, m = m)
    }
    v <- var(unlist(estimators))
    c <- mean(2 * taus + pmax(1, m + 1 - taus))
    ineff.df <- rbind(ineff.df, data.frame(k = k, mfactor = paste0(mfactor[im], " times k"), c = c, v = v))
  }
}

ineff.df$inef <- (ineff.df$c * ineff.df$v) / mcmcvar

library(xtable)
formatted.df <- ineff.df
formatted.df$k <- paste0(formatted.df$k)
cap <- "Cost, variance and inefficiency compared to MCMC, for various choices of $k$ and $m$,
in the baseball batting averages example. \\label{table:baseball}"
colnames(formatted.df) <- c("k", "m", "expected cost", "variance", "inefficiency / MCMC")
formatted.df <- xtable(formatted.df, digits = 4, caption = cap)
formatted.df
print.xtable(formatted.df, include.rownames = FALSE, include.colnames = TRUE, file = "baseball.inefficiency.tex")


# Histogram approximating marginal distributions
k <- 3
m <- 30

nclass <- 40
hist_breaks <- seq(0.1, 0.7, length.out = nclass)
ch3 <- chain[1000:niterations,3]
hist_mcmc <- hist(ch3[ch3>0.1 & ch3 < 0.7], breaks = hist_breaks, plot = F)
hist_mcmc$density
hist_mcmc$mids
# histogram
histogram1 <- histogram_c_chains(c_chains_, 3, k, m, breaks = hist_breaks)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[1])) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red", alpha = 0.5)
g1
ggsave(filename = "baseball.histogram.theta1.pdf", plot = g1, width = 5, height = 5)
#

hist_breaks <- seq(-0.3, 0.85, length.out = nclass)
ch1 <- chain[1000:niterations,1]
hist_mcmc <- hist(ch1[ch1>-0.3 & ch1 < 0.85], breaks = hist_breaks, plot = F)
# hist_mcmc <- hist(chain[1000:niterations,1], nclass = nclass, plot = F)
hist_mcmc$density
hist_mcmc$mids
histogram1 <- histogram_c_chains(c_chains_, 1, k, m, breaks = hist_mcmc$breaks)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(mu)) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red", alpha = 0.5)
g1
ggsave(filename = "baseball.histogram.mu.pdf", plot = g1, width = 5, height = 5)


hist_breaks <- seq(0, 1.2, length.out = nclass)
ch1 <- chain[1000:niterations,2]
hist_mcmc <- hist(ch1[ch1>0 & ch1 < 1.2], breaks = hist_breaks, plot = F)
# hist_mcmc <- hist(chain[1000:niterations,1], nclass = nclass, plot = F)
hist_mcmc$density
hist_mcmc$mids
histogram1 <- histogram_c_chains(c_chains_, 2, k, m, breaks = hist_mcmc$breaks)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(A)) + ylab("density")
g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red", alpha = 0.5)
g1
ggsave(filename = "baseball.histogram.A.pdf", plot = g1, width = 5, height = 5)


