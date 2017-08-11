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
X <- scale(diabetes$x)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

# pb <- get_blasso(Y, X, 0.1)
# niterations <- 5e3
# chain_start <- pb$rinit()
# chain <- matrix(nrow = niterations, ncol = length(chain_start))
# chain[1,] <- chain_start
# for (iteration in 2:niterations){
#   chain[iteration,] <- pb$gibbs_kernel(chain[iteration-1,])
# }
# par(mfrow = c(1,1))
# matplot(chain[,1:10], type = "l")
# c_chain <- coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, max_iterations = 100)
# c_chain$meetingtime
##


lambdas <- 10^(seq(from = -2, to = 3, length.out = 25))
nsamples <- 100
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
  k <- floor(as.numeric(quantile(meetingtimes, probs = 0.90)) + 1)
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
  save(df, nsamples, lambdas, file = "diabetes.cchains.RData")
}
load(file = "diabetes.cchains.RData")


head(df)

# # add another lambda
# lambda <- (lambdas[length(lambdas)-1] - lambdas[length(lambdas)-2]) / 2
# ilambda <- length(lambdas) + 1
# c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#   pb <- get_blasso(Y, X, lambda)
#   coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit)
# }
# meetingtimes <- sapply(c_chains_, function(x) x$meetingtime)
# print(summary(meetingtimes))
# k <- floor(as.numeric(quantile(meetingtimes, probs = 0.90)) + 1)
# K <- floor(10*k)
# cat("found k =", k, "and K =", K, "\n")
# c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
#   pb <- get_blasso(Y, X, lambda)
#   continue_coupled_chains(c_chains_[[irep]], pb$gibbs_kernel, K)
# }
# mean_estimators <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
#   H_bar(c_chains_continued_[[irep]], h = function(x) x[1:p], k = k, K = K)
# }
# sum_estimators <- colSums(mean_estimators)
# sum_squareestimators <- colSums(mean_estimators^2)
# df <- rbind(df, data.frame(ilambda = rep(ilambda, p), lambda = rep(lambda, p), component = 1:p, k = rep(k, p),
#                            K = rep(K, p), nsamples = rep(nsamples, p),
#                            sum_est = matrix(sum_estimators, nrow = p),
#                            sum_squareest = matrix(sum_squareestimators, nrow = p)))
# lambdas <- c(lambdas, lambda)
# save(df, nsamples, lambdas, file = "diabetes.cchains_augmented.RData")


# df$signal <- df$component <= 5
# df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
# df$mean <- (df$sum_est/df$nsamples)
# g <- ggplot(df %>% filter(component <= 10), aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- g + geom_segment(aes(y = mean - 5*sd / sqrt(nsamples), yend = mean + 5*sd / sqrt(nsamples), xend = lambda))
# # g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
# # g <- g + scale_color_colorblind()
# g + scale_x_log10()
# #
# ggplot(df, aes(x = lambda, y = k)) + geom_line() + geom_point() + scale_x_log10()
# ggplot(df, aes(x = lambda, y = sd/abs(mean), group = component, colour = factor(component))) + geom_line() + geom_point() + scale_x_log10()
#
# # g <- g + geom_errorbar(aes(ymin = mean - 2*sd / sqrt(nsamples), ymax = mean + 2*sd / sqrt(nsamples)))
