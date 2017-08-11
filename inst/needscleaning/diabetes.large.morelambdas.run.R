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


load(file = "diabetes.large.cchains.RData")
load(file = "diabetes.large.more.cchains.RData")

head(df)
tail(df)
maxlambda <- max(df$ilambda)
lambdas <- lambdas[1:maxlambda]
nlambdas <- length(lambdas)
new_lambdas <- 10^(seq(from = 3.5, to = 4, length.out = 3))
lambdas <- c(lambdas, new_lambdas)
nsamples <- 100 # detectCores() * 3
for (ilambda in (maxlambda+1):length(lambdas)){
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
  save(df, nsamples, lambdas, file = "diabetes.large.more.cchains.RData")
}
load(file = "diabetes.large.more.cchains.RData")


df$sd <- sqrt(df$sum_squareest/df$nsamples - (df$sum_est/df$nsamples)^2)
df$mean <- (df$sum_est/df$nsamples)
# df$signal <- df$component <= 5
g <- ggplot(df %>% filter(component <= 10), aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
# g <- ggplot(df, aes(x = lambda, y = mean, group = component, colour = factor(component))) + geom_line() + geom_point()
g <- g + geom_segment(aes(y = mean - 2*sd / sqrt(nsamples), yend = mean + 2*sd / sqrt(nsamples), xend = lambda))
# g <- g + scale_color_colorblind()
g <- g + scale_x_log10() + theme(legend.position = "none")
g
