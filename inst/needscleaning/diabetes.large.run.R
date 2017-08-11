# load packages
library(debiasedmcmc)
library(ggthemes)
library(lars)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
# registerDoParallel(cores = 1)
setwd("~/Dropbox/PolyaGammaResults/reproduce/")

data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

#
# C <- matrix(XtX, nrow = p)
# summary(as.numeric(abs(C - XtX)))
# C_inv <- solve(C)
# C_inv_chol <- chol(solve(C))
# C_inv_chol2 <- solve(chol(C))
# C[1:5,1:5]
# C_inv_chol[1:5,1:5]
# C_inv_chol2[1:5,1:5]
# x <- matrix(rnorm(p), nrow = 1)
# sigma2 <- 0.56
# fast_dmvnorm(x, rep(0, p), sigma2 * C)
# fast_dmvnorm_chol_inverse(x, rep(0, p), C_inv_chol2 / sqrt(sigma2))
# chol(solve(C))[1:5,1:5]



#
# pb <- get_blasso(Y, X, 1e-5)
# niterations <- 1e3
# chain_start <- pb$rinit()
# chain <- matrix(nrow = niterations, ncol = length(chain_start))
# chain[1,] <- chain_start
# for (iteration in 2:niterations){
#   chain[iteration,] <- pb$gibbs_kernel(chain[iteration-1,])
# }

# par(mfrow = c(1,1))
# matplot(chain[,1:10], type = "l")
# matplot(chain[,11:20], type = "l")
# matplot(chain[,21:30], type = "l")
# matplot(chain[,31:40], type = "l")
# matplot(chain[,41:50], type = "l")

# acf(chain[,p+1])
# par(mfrow = c(2,2))
# hist(chain[,1])
# hist(chain[,2])
# hist(chain[1000:niterations,p+1])
# hist(chain[1000:niterations,p+2])
##

# pb <- get_blasso(Y, X, 1e-5)
# c_chain <- coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, max_iterations = 1e4, preallocate = 1e2)
# print(c_chain$meetingtime)

# for (irep in 1:20){
#   pb <- get_blasso(Y, X, 1e-5)
#   c_chain <- coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, max_iterations = 1e4)
#   print(c_chain$meetingtime)
# }

lambdas <- 10^(seq(from = -4, to = 2, length.out = 20))
nsamples <- 100 # detectCores() * 3
df <- data.frame()
for (ilambda in 1:length(lambdas)){
  lambda <- lambdas[ilambda]
  print(lambda)
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    pb <- get_blasso(Y, X, lambda)
    coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit)
  }
  meetingtimes <- sapply(c_chains_, function(x) x$meetingtime)
  print(summary(meetingtimes))
  meantau <- mean(meetingtimes)
  mediantau <- median(meetingtimes)
  k <- floor(as.numeric(quantile(meetingtimes, probs = 0.95)) + 1)
  K <- floor(10*k)
  cat("found k =", k, "and K =", K, "\n")
  c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
    pb <- get_blasso(Y, X, lambda)
    continue_coupled_chains(c_chains_[[irep]], pb$gibbs_kernel, K)
  }
  mean_estimators <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    H_bar(c_chains_continued_[[irep]], h = function(x) x[1:p], k = k, K = K)
  }
  sum_estimators <- colSums(mean_estimators)
  sum_squareestimators <- colSums(mean_estimators^2)
  df <- rbind(df, data.frame(ilambda = rep(ilambda, p), lambda = rep(lambda, p), component = 1:p, k = rep(k, p),
                             K = rep(K, p), nsamples = rep(nsamples, p),
                             meantau = meantau,
                             mediantau = mediantau,
                             sum_est = matrix(sum_estimators, nrow = p),
                             sum_squareest = matrix(sum_squareestimators, nrow = p)))
  save(df, nsamples, lambdas, file = "diabetes.large.cchains.RData")
}
load(file = "diabetes.large.cchains.RData")


head(df)
tail(df)

df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
df$mean <- (df$sum_est/df$nsamples)
# df$signal <- df$component <= 5
g <- ggplot(df %>% filter(component <= 10), aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- ggplot(df, aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
# g <- g + scale_color_colorblind()
g <- g + scale_x_log10() + theme(legend.position = "none")
g
