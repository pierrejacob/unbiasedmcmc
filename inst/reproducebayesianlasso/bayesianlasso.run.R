# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())

data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)
lambdas <- 10^(seq(from = -2, to = 3, length.out = 25))

# lambda <- lambdas[5]
# pb <- get_blasso(Y, X, lambda)
# unbiasedestimator(single_kernel, coupled_kernel, pb$rinit, h = function(x) x, k = 0, m = 0)

nrep <- 100
df <- data.frame()
for (ilambda in 1:length(lambdas)){
  print(ilambda)
  lambda <- lambdas[ilambda]
  pb <- get_blasso(Y, X, lambda)
  single_kernel <- function(state) list(state = pb$gibbs_kernel(state))
  coupled_kernel <- function(state1, state2){
    res <- pb$coupled_gibbs_kernel(state1, state2)
    return(list(state1 = res$chain_state1,
                state2 = res$chain_state2))
  }
  ues <- foreach(irep = 1:nrep) %dopar% {
     unbiasedestimator(single_kernel, coupled_kernel, pb$rinit, h = function(x) x, k = 0, m = 0)
  }
  meetingtimes <- sapply(ues, function(x) x$meetingtime)
  df <- rbind(df, data.frame(ilambda = ilambda, lambda = lambda, rep = 1:nrep, meetingtime = meetingtimes))
  save(df, lambdas, nrep, file = "bayesianlasso.meetings.RData")
}
save(df, lambdas, nrep, file = "bayesianlasso.meetings.RData")



load("bayesianlasso.meetings.RData")
head(df)
meetingquantile <- df %>% group_by(ilambda) %>% summarise(mean = mean(meetingtime), q99 = quantile(meetingtime, probs = 0.99))
meetingquantile$q99


df <- data.frame()
for (ilambda in 1:length(lambdas)){
  print(ilambda)
  lambda <- lambdas[ilambda]
  k <- meetingquantile$q99[ilambda]
  m <- 10 * k
  pb <- get_blasso(Y, X, lambda)
  single_kernel <- function(state) list(state = pb$gibbs_kernel(state))
  coupled_kernel <- function(state1, state2){
    res <- pb$coupled_gibbs_kernel(state1, state2)
    return(list(state1 = res$chain_state1,
                state2 = res$chain_state2))
  }
  ues <- foreach(irep = 1:nrep) %dopar% {
    unbiasedestimator(single_kernel, coupled_kernel, pb$rinit, h = function(x) x, k = k, m = m)
  }
  # meetingtimes <- sapply(ues, function(x) x$meetingtime)
  for (irep in 1:nrep){
    df <- rbind(df, data.frame(ilambda = ilambda, lambda = lambda, rep = rep(irep, (2*p+1)),
                               component = 1:(2*p+1), estimator = ues[[irep]]$uestimator))
  }
  save(df, lambdas, nrep, file = "bayesianlasso.estimators.RData")
}
save(df, lambdas, nrep, file = "bayesianlasso.estimators.RData")

load("bayesianlasso.estimators.RData")
summary.df <- df %>% group_by(ilambda, lambda, component) %>% summarise(m = mean(estimator), sd = sd(estimator), nsamples = n())
summary.df$component %>% unique
summary.df$nsamples %>% unique

g <- ggplot(summary.df %>% filter(component <= 64), aes(x = lambda, y = m, group = component)) + geom_line() + scale_x_log10()
g + geom_segment(aes(y = m - 2*sd / sqrt(nrep), yend = m + 2*sd / sqrt(nrep), xend = lambda))

# g <- ggplot(summary.df %>% filter(component <= 64), aes(x = ilambda, y = m, group = component)) + geom_line() + scale_x_log10()
# g + geom_segment(aes(y = m - 2*sd / sqrt(nsamples), yend = m + 2*sd / sqrt(nsamples), xend = ilambda))


# compute extra estimators for first 10 lambdas
nrep <- 1000
extra.df <- data.frame()
for (ilambda in 1:10){
  print(ilambda)
  lambda <- lambdas[ilambda]
  k <- meetingquantile$q99[ilambda]
  m <- 10 * k
  pb <- get_blasso(Y, X, lambda)
  single_kernel <- function(state) list(state = pb$gibbs_kernel(state))
  coupled_kernel <- function(state1, state2){
    res <- pb$coupled_gibbs_kernel(state1, state2)
    return(list(state1 = res$chain_state1,
                state2 = res$chain_state2))
  }
  ues <- foreach(irep = 1:nrep) %dopar% {
    unbiasedestimator(single_kernel, coupled_kernel, pb$rinit, h = function(x) x, k = k, m = m)
  }
  # meetingtimes <- sapply(ues, function(x) x$meetingtime)
  for (irep in 1:nrep){
    extra.df <- rbind(extra.df, data.frame(ilambda = ilambda, lambda = lambda, rep = rep(irep, (2*p+1)),
                               component = 1:(2*p+1), estimator = ues[[irep]]$uestimator))
  }
  save(extra.df, file = "bayesianlasso.extraestimators.RData")
}
save(extra.df, file = "bayesianlasso.extraestimators.RData")

load("bayesianlasso.extraestimators.RData")

summary.extra.df <- rbind(df, extra.df) %>% group_by(ilambda, lambda, component) %>% summarise(m = mean(estimator), sd = sd(estimator), nsamples = n())
summary.extra.df$nsamples %>% unique

g <- ggplot(summary.extra.df %>% filter(component <= 64), aes(x = lambda, y = m, group = component)) + geom_line() + scale_x_log10()
g + geom_segment(aes(y = m - 2*sd / sqrt(nsamples), yend = m + 2*sd / sqrt(nsamples), xend = lambda))

# g <- ggplot(summary.extra.df %>% filter(component <= 64), aes(x = ilambda, y = m, group = component)) + geom_line() + scale_x_log10()
# g + geom_segment(aes(y = m - 2*sd / sqrt(nsamples), yend = m + 2*sd / sqrt(nsamples), xend = ilambda))
