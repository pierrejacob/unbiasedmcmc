# This script plays with the computation of estimators
# using signed measure representations

library(unbiasedmcmc)
library(doParallel)
library(doRNG)
library(dplyr)
registerDoParallel(cores = detectCores()-2)
rm(list = ls())
set.seed(1)

#### Bivariate target
target <- function(x) fast_dmvnorm(matrix(x, nrow = 1), mean = c(0,0.5), covariance = diag(c(0.5, 0.2)))
Sigma_proposal <- diag(1, 2)
rinit <- function(){
  x <- fast_rmvnorm(1, c(10, 5), diag(3, 2, 2))
  return(list(chain_state = x, current_pdf = target(x)))
}
kernels <- get_mh_kernels(target, Sigma_proposal)

cx <- sample_coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit)
cx$samples1 %>% head
cx$meetingtime
dim(cx$samples1)[1]
dim(cx$samples1)[2]

ms <- c_chains_to_measure_as_list(cx, 0, 0)
tail(ms$atoms)
cat(sum(ms$weights * ms$atoms[,1]), sum(ms$weights * ms$atoms[,ncol(ms$atoms)]), "\n")
H_bar(cx, h = function(x) x, k = 0, m = 0)

ms$MCMC
sum(ms$weights[ms$MCMC==1])
sum(ms$weights[ms$MCMC==0])

cx <- sample_coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit, m = 100)
cx$meetingtime
ms <- c_chains_to_measure_as_list(cx, 30, 100)
tail(ms$atoms)
cat(sum(ms$weights * ms$atoms[,1]), sum(ms$weights * ms$atoms[,ncol(ms$atoms)]), "\n")
H_bar(cx, h = function(x) x, k = 30, m = 100)
ms$MCMC
sum(ms$weights[ms$MCMC])
sum(ms$weights[!ms$MCMC])

cx <- sample_coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit, m = 100, lag = 5)
cx$meetingtime
ms <- c_chains_to_measure_as_list(cx, 3, 100)
ms$MCMC
sum(ms$weights[ms$MCMC])
sum(ms$weights[!ms$MCMC])

cat(sum(ms$weights * ms$atoms[,1]), sum(ms$weights * ms$atoms[,ncol(ms$atoms)]), "\n")
H_bar(cx, h = function(x) x, k = 3, m = 100)

c_chains_to_dataframe(cx, k = 3, m = 100)

nrepeats <- 10
coupledchains_ <- foreach(rep = 1:nrepeats) %dorng% {
  sample_coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit, m = 100, lag = 1)
}
summary(sapply(coupledchains_, function(l) l$meetingtime))

k <- 20
m <- 100
df_ <- c_chains_to_dataframe(coupledchains_, k, m, dopar = T)
head(df_)
colSums(df_$weight * df_[,4:ncol(df_)])
rowMeans(sapply(X = coupledchains_, FUN = function(x) H_bar(x, h = function(v) v, k = k, m = m)))

cat(sum(df_$weight * df_$atom.1 * cos(df_$atom.2)), "\n")
mean(sapply(X = coupledchains_, FUN = function(x) H_bar(x, h = function(v) v[1] * cos(v[2]), k = k, m = m)))




##### Univariate target
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-2, 2), sd = 0.5, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}
Sigma_proposal <- diag(2, 1, 1)
rinit <- function(){
  x <- rnorm(1, 5, 5)
  return(list(chain_state = x, current_pdf = target(x)))
}
kernels <- get_mh_kernels(target, Sigma_proposal)

nrepeats <- 10000
m <- 200
coupledchains_ <- foreach(rep = 1:nrepeats) %dorng% {
  sample_coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit, m = m, lag = 50)
}
summary(sapply(coupledchains_, function(l) l$meetingtime))
k <- 20
df_ <- c_chains_to_dataframe(coupledchains_, k, m, dopar = T)
head(df_)

sum(df_$weight * df_$atom.1)
colSums(df_$weight * df_[,4:ncol(df_),drop=F])
mean(sapply(X = coupledchains_, FUN = function(x) H_bar(x, h = function(v) v, k = k, m = m)))

df_ %>% summarise(mean1 = sum(weight * atom.1))
df_ %>% group_by(MCMC) %>% summarise(mean1 = sum(weight * atom.1), sumweight = sum(weight)) %>% ungroup() %>% as.data.frame
mean(df_$MCMC==1)
MCMC_traj <- sapply(coupledchains_, function(x) x$samples1[k:m,])
mean(MCMC_traj)


hist1 <- histogram_c_chains(coupledchains_, component = 1, k = 150, m = m)
setmytheme()
plot_histogram(hist1)



dim(MCMC_traj)

matplot(t(MCMC_traj[,1:100]), type = "l")
