library(debiasedmcmc)
rm(list = ls())
# double well potential log-pdf
target <- function(x) -0.1 * ( ((x[1]-1)^2-x[2]^2 )^2 + 10*(x[1]^2-5)^2+ (x[1]+x[2])^4 + (x[1]-x[2])^4)
# compute marginals pdfs and nornalizing constant
targetmarginal1 <- function(x) integrate(function(a) sapply(a, function(z) exp(target(c(x, z)))), lower = -6, upper = 6, subdivisions = 1e4)$value
targetmarginal2 <- function(x) integrate(function(a) sapply(a, function(z) exp(target(c(z, x)))), lower = -6, upper = 6, subdivisions = 1e4)$value
Z <- integrate(function(x) sapply(x, function(z) targetmarginal1(z)), lower = -6, upper = 6, subdivisions = 1e4)$value

# xgrid <- seq(from = -4, to = 4, length.out = 1e2)
# ygrid <- seq(from = -4, to = 4, length.out = 1e2)
# df_ <- expand.grid(xgrid, ygrid)
# df_$z <- apply(df_, 1, target)
# library(ggplot2)
# ggplot(df_, aes(x = Var1, y = Var2, z = z)) + geom_contour(binwidth = 5)

# define random walk Metropolis-Hastings kernels
kernels <- get_mh_kernels(target, Sigma_proposal = diag(1, 2, 2))
# define initial distribution
rinit <- function(){
  x <- rnorm(n = 2, mean = 2, sd = 25)
  return(list(chain_state = x, current_pdf = target(x)))
}
## register parallel cores
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
nrepeats <- 5e2
# sample meeting times "tau"
meetings_ <- foreach(rep = 1:nrepeats) %dorng% {
  sample_meetingtime(kernels$single_kernel, kernels$coupled_kernel, rinit)
}
hist(sapply(meetings_, function(x) x$meetingtime), nclass = 50, xlab = "meeting time", main = "")

## now run coupled chains for max(m, tau) steps
coupledchains_ <- foreach(rep = 1:nrepeats) %dorng% {
  sample_coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit, m = 10000, lag = 1000)
}
names(coupledchains_[[1]])
hist(sapply(coupledchains_, function(x) x$meetingtime))

# let's see what happened to the chains that took longest to meet
index_max <- which.max(sapply(coupledchains_, function(x) x$meetingtime))
plot(coupledchains_[[index_max]]$samples1[,1], coupledchains_[[index_max]]$samples1[,2], type = "l", xlim = c(-5, 5), ylim = c(-5, 5))
lines(coupledchains_[[index_max]]$samples2[,1], coupledchains_[[index_max]]$samples2[,2], col = "red")
matplot(coupledchains_[[index_max]]$samples1, type = "l", col = "black")
matplot(coupledchains_[[index_max]]$samples2, type = "l", col = "red", add = T)

library(ggplot2)
hist1 <- histogram_c_chains(coupledchains_, component = 1, k = 2000, m = 5000, dopar = FALSE)
plot_histogram(hist1) + stat_function(fun = function(x) sapply(x, function(v) targetmarginal1(v)/Z)) + xlim(-4,4)

hist2 <- histogram_c_chains(coupledchains_, component = 2, k = 2000, m = 10000, dopar = TRUE)
plot_histogram(hist2) + stat_function(fun = function(x) sapply(x, function(v) targetmarginal2(v)/Z)) + xlim(-4,4)


