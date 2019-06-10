# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
library(dplyr)

registerDoParallel(cores = detectCores()-2)

data(diabetes)
X <- scale(diabetes$x2)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

## Modify lambdas
# lambdas <- 10^(seq(from = -2, to = 3, length.out = 25))
lambdas <- 10^(seq(from = -2, to = 1, length.out = 15))

# unbiasedestimator <- function(single_kernel, coupled_kernel, rinit, ..., h = function(x) x, k = 0, m = 1, max_iterations = Inf){
#   chain_state1 <- rinit()
#   chain_state2 <- rinit()
#   # mcmcestimator computes the sum of h(X_t) for t=k,...,m
#   mcmcestimator <- h(chain_state1)
#   dimh <- length(mcmcestimator)
#   if (k > 0){
#     mcmcestimator <- rep(0, dimh)
#   }
#   # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
#   correction <- rep(0, dimh)
#   chain_state1 <- single_kernel(chain_state1, ...)$state
#   # chain_state1 <- sres1$state
#   if (k == 0){
#     correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
#   }
#   if (k <= 1 && m >= 1){
#     mcmcestimator <- mcmcestimator + h(chain_state1)
#   }
#   # current_nsamples1 <- current_nsamples1 + 1
#   # samples1[current_nsamples1,] <- chain_state1
#   iter <- 1
#   meet <- FALSE
#   finished <- FALSE
#   meetingtime <- Inf
#   # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
#   while (!finished && iter < max_iterations){
#     iter <- iter + 1
#     if (meet){
#       chain_state1 <- single_kernel(chain_state1, ...)$state
#       # chain_state1 <- sres1$state
#       chain_state2 <- chain_state1
#       if (k <= iter && iter <= m){
#         mcmcestimator <- mcmcestimator + h(chain_state1)
#       }
#     } else {
#       res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, ...)
#       chain_state1 <- res_coupled_kernel$state1
#       chain_state2 <- res_coupled_kernel$state2
#       if (all(chain_state1 == chain_state2) && !meet){
#         # recording meeting time tau
#         meet <- TRUE
#         meetingtime <- iter
#       }
#       if (k <= iter){
#         if (iter <= m){
#           mcmcestimator <- mcmcestimator + h(chain_state1)
#         }
#         correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
#       }
#     }
#     # stop after max(m, tau) steps
#     if (iter >= max(meetingtime, m)){
#       finished <- TRUE
#     }
#   }
#   mcmcestimator <- mcmcestimator / (m - k + 1)
#   uestimator <- mcmcestimator + correction
#   return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
#               meetingtime = meetingtime, iteration = iter, finished = finished))
# }


# lambda <- lambdas[5]
# pb <- get_blasso(Y, X, lambda)
# single_kernel <- function(state) list(state = pb$gibbs_kernel(state))
# coupled_kernel <- function(state1, state2){
#   res <- pb$coupled_gibbs_kernel(state1, state2)
#   return(list(state1 = res$chain_state1,
#               state2 = res$chain_state2))
# }
# unbiasedestimator(single_kernel, coupled_kernel, pb$rinit, h = function(x) x, k = 0, m = 0)

nrep <- 100
df <- data.frame()
for (ilambda in 1:length(lambdas)){
  print(ilambda)
  lambda <- lambdas[ilambda]
  pb <- get_blasso(Y, X, lambda)
  # single_kernel <- function(state) list(state = pb$gibbs_kernel(state))
  # coupled_kernel <- function(state1, state2){
  #   res <- pb$coupled_gibbs_kernel(state1, state2)
  #   return(list(state1 = res$chain_state1,
  #               state2 = res$chain_state2))
  # }
  meetings <- foreach(irep = 1:nrep) %dorng% {
    sample_meetingtime(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit)
  }
  meetingtimes <- sapply(meetings, function(x) x$meetingtime)
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
  # single_kernel <- function(state) list(state = pb$gibbs_kernel(state))
  # coupled_kernel <- function(state1, state2){
  #   res <- pb$coupled_gibbs_kernel(state1, state2)
  #   return(list(state1 = res$chain_state1,
  #               state2 = res$chain_state2))
  # }
  ues <- foreach(irep = 1:nrep) %dorng% {
    sample_unbiasedestimator(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, h = function(x) x, k = k, m = m)
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
  ues <- foreach(irep = 1:nrep) %dorng% {
    sample_unbiasedestimator(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, h = function(x) x, k = k, m = m)
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
