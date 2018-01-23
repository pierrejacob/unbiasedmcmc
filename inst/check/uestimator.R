library(debiasedmcmc)
setwd("~/Dropbox/UnbiasedMCMCResults/reproduce/")
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores() - 1)
#
# simulate data
n <- 500
p <- 1000
SNR <- 3
s_star <- 10
s0 <- 100
sigma0 <- 1
beta_star <- SNR * sqrt(sigma0^2 * log(p) / n) * c(2,-3,2,2,-3,3,-2,3,-2,3, rep(0, p-10))
# independent design
X <- matrix(rnorm(n * p), nrow = n, ncol = p) # fast_rmvnorm_chol(n, rep(0, p), diag(1, p, p))
X <- scale(X)
Y <- X %*% matrix(beta_star, ncol = 1) + rnorm(n, 0, sigma0)
Y <- scale(Y)
Y2 <- (t(Y) %*% Y)[1,1]
g <- p^3
kappa <- 1
proportion_singleflip <- 0.5
vs <- get_variableselection(Y,X,g,kappa,s0,proportion_singleflip)
rinit <- vs$rinit
single_kernel <- vs$single_kernel
coupled_kernel <- vs$coupled_kernel
coupled_chains <- vs$coupled_chains
marginal_likelihood <- vs$marginal_likelihood
prior <- vs$prior
#



k <- 0
m <- 0
h <- function(x) x
set.seed(17)
cc <- coupled_chains(single_kernel, coupled_kernel, rinit, m = m)
cc$meetingtime
h_bar_estimator <- H_bar(cc, h = h, k = k, m = m)
#
set.seed(17)
ue <- unbiasedestimator_vs(single_kernel, coupled_kernel, rinit, h = h, k = k, m = m)
ue$meetingtime
# comparison
summary(ue$uestimator - h_bar_estimator)

