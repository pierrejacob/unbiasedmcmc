# load packages
library(debiasedmcmc)
library(ggthemes)
library(lars)
setmytheme()
rm(list = ls())
set.seed(28)
registerDoParallel(cores = detectCores())
# registerDoParallel(cores = 1)
setwd("~/Dropbox/PolyaGammaResults/reproduce/")

data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)


# load(file = "diabetes.large.cchains.RData")
load(file = "diabetes.large.more.cchains.RData")

head(df)
tail(df)
lambdas

unique(df$lambda)
unique(df$ilambda)
which(unique(df$lambda) < 1e-1)

# so let's add more samples to first 10 values

nsamples <- 1000
for (ilambda in 1:10){
  print(ilambda)
  lambda <- lambdas[ilambda]
  df_sub <- df[df$ilambda==ilambda,]
  kK <- (df_sub %>% select(k,K) %>% colMeans)
  k <- as.numeric(kK[1])
  K <- as.numeric(kK[2])
  old_nsamples <- as.numeric(df_sub %>% select(nsamples) %>% colMeans)
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    pb <- get_blasso(Y, X, lambda)
    coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, K = K)
  }

  meetingtimes <- sapply(c_chains_, function(x) x$meetingtime)
  print(summary(meetingtimes))
  mean_estimators <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    H_bar(c_chains_[[irep]], h = function(x) x[1:p], k = k, K = K)
  }

  sum_estimators <- colSums(mean_estimators)
  sum_squareestimators <- colSums(mean_estimators^2)

  df_sub$sum_est <- df_sub$sum_est + sum_estimators
  df_sub$sum_squareest <- df_sub$sum_squareest + sum_squareestimators
  df_sub$nsamples <- df_sub$nsamples + nsamples

  df[df$ilambda==ilambda,] <-  df_sub
  save(df, lambdas, file = "diabetes.large.more.more.cchains.RData")
}
load(file = "diabetes.large.more.more.cchains.RData")


# df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
# df$mean <- (df$sum_est/df$nsamples)
# # df$signal <- df$component <= 5
# g <- ggplot(df %>% filter(component <= 10), aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# # g <- ggplot(df, aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
# # g <- g + scale_color_colorblind()
# g <- g + scale_x_log10() + theme(legend.position = "none")
# g
