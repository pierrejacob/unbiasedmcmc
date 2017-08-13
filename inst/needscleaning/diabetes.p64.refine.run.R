# load packages
library(debiasedmcmc)
library(ggthemes)
# library(lars)
library(parallel)
setmytheme()
rm(list = ls())
set.seed(28)
# registerDoParallel(cores = detectCores())
setwd("/home/pierre/Dropbox/PolyaGammaResults")
data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

load(file = "diabetes.cchains.p64.RData")
lambdas

unique(df$lambda)
unique(df$ilambda)
which(unique(df$lambda) < 1e-1)

# so let's add more samples to first 10 values
# cl <- makeCluster(detectCores(), type = "FORK")
# clusterSetRNGStream(cl, 12)
# clusterEvalQ(cl, library("debiasedmcmc"))
# clusterEvalQ(cl, library("lars"))
# clusterEvalQ(cl, {data(diabetes);
#   X <- scale(diabetes$x2);
#   Y <- matrix(scale(diabetes$y), ncol = 1);
# })
f <- function(x, lambda, K = 1){
  pb <- get_blasso(Y, X, lambda);
  return(coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, K = K))
}
# clusterExport(cl, "f")
#

nsamples <- 20
for (ilambda in 1:10){
  print(ilambda)
  lambda <- lambdas[ilambda]
  # clusterExport(cl, "lambda")
  df_sub <- df[df$ilambda==ilambda,]
  kK <- (df_sub %>% select(k,K) %>% colMeans)
  k <- as.numeric(kK[1])
  K <- as.numeric(kK[2])
  old_nsamples <- as.numeric(df_sub %>% select(nsamples) %>% colMeans)
  # clusterExport(cl, "K")
  # c_chains_ <- parLapply(cl, 1:nsamples, function(x) f(x, lambda, K))
  c_chains_ <- list()
  mean_estimators <- matrix(nrow = nsamples, ncol = p)
  for (isample in 1:nsamples){
    c_chains_[[isample]] <- f(1, lambda, K)
    mean_estimators[isample,] <- H_bar(c_chains_[[isample]], h = function(x) x[1:p], k = k, K = K)
  }
  meetingtimes <- sapply(c_chains_, function(x) x$meetingtime)
  print(summary(meetingtimes))
  # mean_estimators <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
  #   H_bar(c_chains_[[irep]], h = function(x) x[1:p], k = k, K = K)
  # }
  sum_estimators <- colSums(mean_estimators)
  sum_squareestimators <- colSums(mean_estimators^2)
  df_sub$sum_est <- df_sub$sum_est + sum_estimators
  df_sub$sum_squareest <- df_sub$sum_squareest + sum_squareestimators
  df_sub$nsamples <- df_sub$nsamples + nsamples
  df[df$ilambda==ilambda,] <-  df_sub
  save(df, lambdas, file = "diabetes.cchains.p64.refine.RData")
}
load(file = "diabetes.cchains.p64.refine.RData")


# df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
# df$mean <- (df$sum_est/df$nsamples)
# # df$signal <- df$component <= 5
# g <- ggplot(df %>% filter(component <= 10), aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# # g <- ggplot(df, aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
# # g <- g + scale_color_colorblind()
# g <- g + scale_x_log10() + theme(legend.position = "none")
# g
