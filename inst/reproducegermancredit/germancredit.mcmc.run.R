
# load packages
library(debiasedmcmc)
rm(list = ls())
set.seed(21)

#
##  This example is about the Polya Gamma Gibbs sampler for logistic regression models, as applied to the German credit data of Lichman 2013.

data(germancredit)
n <- nrow(X)
p <- ncol(X)

# prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
logistic_setting <- logisticregression_precomputation(Y, X, b, B)
# define MCMC transition kernel
single_kernel <- function(chain_state, logistic_setting){
  zs <- abs(logisticregression_xbeta(logistic_setting$X, t(chain_state)))
  w <- BayesLogit::rpg(logistic_setting$n, h=1, z=zs)
  res <- logisticregression_m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(chain_state)
}

rinit <- function(){
  t(fast_rmvnorm(1, mean = b, covariance = B))
}

chain_state <- rinit()
zs <- abs(logisticregression_xbeta(logistic_setting$X, t(chain_state)))
w <- BayesLogit::rpg(logistic_setting$n, h=1, z=zs)
res <- logisticregression_m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
Sigma_ <- solve(res$Sigma_inverse)
Sigma_alt <- t(res$Cholesky) %*% res$Cholesky

Sigma_[1:5,1:5]
Sigma_alt[1:5,1:5]

mean_ <- res$m
n <- 1e5
xx_2 <- fast_rmvnorm_chol(n, mean_, res$Cholesky)
colMeans(xx_2[,1:5])
mean_[1:5]
cov(xx_2[,1:5])
Sigma_[1:5,1:5]
# xx_3 <- fast_rmvnorm_chol(n, mean_, chol(Sigma_))
# colMeans(xx_3[,1:5])
# cov(xx_3[,11:15])
# Sigma_[11:15,11:15]
# modify niterations

niterations <- 10000

chain <- matrix(nrow = niterations, ncol = p)
chain[1,] <- rinit()
for (iteration in 2:niterations){
  chain[iteration,] <- single_kernel(chain[iteration-1,], logistic_setting)
}
save(niterations, chain, file = "germancredit.mcmc.RData")
load("germancredit.mcmc.RData")
matplot(chain[100:niterations,1:5], type = "l")

idx <- which('Instalment.per.cent' == colnames(X))
hist(chain[100:niterations,idx])
