library(debiasedmcmc)
library(doParallel)
library(doRNG)
# register parallel cores
registerDoParallel(cores = detectCores()-2)
setmytheme()
rm(list = ls())
set.seed(21)

target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-2, 2), sd = 0.2, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

kernels <- get_mh_kernels(target, 2^2)
single_kernel <- kernels$single_kernel
coupled_kernel <- kernels$coupled_kernel
rinit <- function(){
  x <- rnorm(1)
  list(chain_state = x, current_pdf = target(x))
}
c_chains <- sample_coupled_chains(single_kernel, coupled_kernel, rinit)


nsamples <- 1000
meetingtime <- foreach(irep = 1:nsamples, .combine = c) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit)$meetingtime
}
summary(meetingtime)
hist(meetingtime)
##
m <- 100
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m, lag = 50)
}
k <- 20
index_ <- which(sapply(c_chains_, function(x) x$meetingtime) > 30)[1]
c_chains_[[index_]]
H_bar(c_chains_[[index_]], h = function(x) x, k = k, m = m)
#
summary(sapply(c_chains_, function(x) x$meetingtime))
# histogram
component <- 1
hist1 <- histogram_c_chains(c_chains_, component, k, m)
barplot(height = hist1$proportions)
plot_histogram(hist1)
library(ggplot2)
df_ <- data.frame(x = hist1$mids, y = hist1$proportions / hist1$width)
g <- ggplot(df_, aes(x = x, xend = x, y = 0, yend = y))
g <- g + geom_segment(aes(x = x, xend = x, y = 0, yend = y)) + ylab("density")
g
