# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
setwd("~/Dropbox/PolyaGammaResults/pgg")
#

# Generate data
n <- 1000 # number of observations
p <- 100 # number of covariates
## create covariates
X <- fast_rmvnorm(n, rep(0, p), diag(1, nrow = p, ncol = p))
# create outcome
beta_star <- 1:p / p # data-generating parameters
Y <- rbinom(n, 1, expit(apply(X = X, MARGIN = 1, FUN = function(row) sum(row * beta_star))))
# specify prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(50, p, p); B[1,1] <- 800
# precompute some quantities useful both for the Polya-Gamma sampler
# and for the unbiased estimators
logistic_setting <- logistic_precomputation(Y, X, b, B)
## Polya-Gamma sampler

single_kernel <- function(chain_state){
  beta <- chain_state[1:p]
  zs <- abs(xbeta(X, beta))
  w <- rpg(n, h=1, z=zs)
  res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  beta <- fast_rmvnorm_chol(1, res$m, res$Cholesky)
  return(beta)
}

niterations <- 10000
betas <- matrix(0, ncol=p, nrow=niterations)
beta <- rep(0, p)
for (i in 1:niterations){
  beta <- single_kernel(beta)
  betas[i,] <- beta
}
# save(betas, file = "pgg.mcmc.RData")
# load("pgg.mcmc.RData")
matplot(betas[,1:10], type = "l")
acf(betas[,1])

sample_beta_dist <- function(w1, w2, logistic_setting){
  X <- logistic_setting$X
  invB <- logistic_setting$invB
  res1 <- m_and_sigma(w1, X, invB, logistic_setting$KTkappaplusinvBtimesb)
  res2 <- m_and_sigma(w2, X, invB, logistic_setting$KTkappaplusinvBtimesb)
  distance <- sqrt((sum(res1$m - res2$m)^2))
  cat("distance:", distance, "\n")
  if(distance <= 10){
    x <- gaussian_max_coupling_cholesky_R(res1$m, res2$m, t(res1$Cholesky), t(res2$Cholesky), t(res1$Cholesky_inverse), t(res2$Cholesky_inverse))
    beta1 <- x[,1]
    beta2 <- x[,2]
  } else {
    x <- gaussian_opt_transport(1, res1$m, res2$m, t(res1$Cholesky), t(res2$Cholesky), t(res1$Cholesky_inverse), t(res2$Cholesky_inverse))
    beta1 <- x[[1]][,1]
    beta2 <- x[[1]][,2]
  }
  return(list(beta1=beta1, beta2=beta2))
}

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2){
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
  betas <- sample_beta_dist(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  return(list(chain_state1=cbind(betas$beta1),
              chain_state2=cbind(betas$beta2)))
}

rinit <- function() fast_rmvnorm(1, rep(0, p), diag(1, p, p))
# rinit <- function() fast_rmvnorm(1, colMeans(betas), cov(betas))
c_chain <- coupled_chains(single_kernel, coupled_kernel, rinit)
c_chain$meetingtime

nsamples <- 50
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit)
}

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)

x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g
# ggsave(filename = "propensity.meetingtime.pdf", plot = g, width = 5, height = 5)

k <- 1000
K <- 2000

# ##
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}
sum(sapply(c_chains_continued_, function(x) x$iteration))

# histogram
k <- 950
K <- 1000
histogram1 <- histogram_c_chains(c_chains_continued_, 33, k, K, nclass = 50)
g1 <- plot_histogram(histogram1, with_bar = F) + xlab(expression(theta)) + ylab("density")
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

summary(abs(est_mean - colMeans(betas))/abs(colMeans(betas)))
