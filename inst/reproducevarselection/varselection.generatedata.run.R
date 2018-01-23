library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(1)
# registerDoParallel(cores = detectCores())

# simulate data
n <- 1000
p <- 5000
SNR <- 1
s_star <- 10
s0 <- 100
sigma0 <- 1
beta_star <- SNR * sqrt(sigma0^2 * log(p) / n) * c(2,-3,2,2,-3,3,-2,3,-2,3, rep(0, p-10))
# independent design
X <- matrix(rnorm(n * p), nrow = n, ncol = p) # fast_rmvnorm_chol(n, rep(0, p), diag(1, p, p))
X <- scale(X)
Y <- X %*% matrix(beta_star, ncol = 1) + rnorm(n, 0, sigma0)
Y <- scale(Y)

save(X, Y, beta_star, file = paste0("varselection.dataSNR", SNR, ".RData"))
