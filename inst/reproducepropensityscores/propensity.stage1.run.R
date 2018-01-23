# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores())
#

# Generate data
n <- 1000 # number of observations
p <- 6 # number of covariates
## create covariates
X <- fast_rmvnorm(n, rep(0, p), diag(1, nrow = p, ncol = p))
# create outcome
beta_star <- 1:p / 10 # data-generating parameters
Y <- rbinom(n, 1, expit(apply(X = X, MARGIN = 1, FUN = function(row) sum(row * beta_star))))
# specify prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(50, p, p); B[1,1] <- 800
# precompute some quantities useful both for the Polya-Gamma sampler
# and for the unbiased estimators
logistic_setting <- logistic_precomputation(Y, X, b, B)
save(logistic_setting, file = "propensity.data.RData")

## Polya-Gamma sampler
single_kernel <- function(chain_state){
  beta <- chain_state[1:p]
  zs <- abs(xbeta(X, beta))
  w <- rpg(n, h=1, z=zs)
  res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  beta <- fast_rmvnorm_chol(1, res$m, res$Cholesky)
  return(beta)
}

# niterations <- 1000
# betas <- matrix(0, ncol=p, nrow=niterations)
# beta <- rep(0, p)
# for (i in 1:niterations){
#   beta <- single_kernel(beta)
#   betas[i,] <- beta
# }
# matplot(betas, type = "l")

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2){
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  return(list(chain_state1=cbind(betas$beta1),
              chain_state2=cbind(betas$beta2)))
}
rinit <- function() fast_rmvnorm(1, rep(0, p), diag(1, p, p))
nsamples <- 10000

filename <- "propensity.stage1.c_chains.RData"
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit)
}
save(nsamples, c_chains_, file = filename)
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)

x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g
# ggsave(filename = "propensity.meetingtime.pdf", plot = g, width = 5, height = 5)


k <- 20
K <- k

# ##
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}
save(k, K, c_chains_, c_chains_continued_, file = filename)

sum(sapply(c_chains_continued_, function(x) x$iteration))

# histogram
histogram1 <- histogram_c_chains(c_chains_continued_, 3, k, K, nclass = 20)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta)) + ylab("density")
g1


mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], k = k, K = K)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], h = function(x) x^2, k = k, K = K)
}

est_mean <- rep(0, p)
est_var <- rep(0, p)

for (component in 1:p){
  cat("component ", component, "\n")
  estimators <- sapply(mean_estimators, function(x) x[component])
  cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
  s_estimators <- sapply(square_estimators, function(x) x[component])
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
  est_mean[component] <- mean(estimators)
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
}

