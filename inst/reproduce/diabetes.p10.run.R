# load packages
library(debiasedmcmc)
library(ggthemes)
library(parallel)

setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
cl <- makeCluster(detectCores(), type = "FORK")
clusterSetRNGStream(cl, 12)
clusterEvalQ(cl, library("debiasedmcmc"))
clusterEvalQ(cl, library("lars"))
clusterEvalQ(cl, {data(diabetes);
  X <- scale(diabetes$x);
  Y <- matrix(scale(diabetes$y), ncol = 1);
})


data(diabetes)
X <- scale(diabetes$x)
Y <- matrix(scale(diabetes$y), ncol = 1)
p <- ncol(X)

f <- function(x, lambda, K = 1){
  pb <- get_blasso(Y, X, lambda);
  return(coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit, K = K))
}
clusterExport(cl, "f")


lambdas <- 10^(seq(from = -2, to = 3, length.out = 25))
nsamples <- 100
df <- data.frame()
for (ilambda in 1:length(lambdas)){
  lambda <- lambdas[ilambda]
  print(lambda)
  clusterExport(cl, "lambda")
  c_chains_ <- parLapply(cl, 1:nsamples, function(x) f(x, lambda))
  meetingtimes <- sapply(c_chains_, function(x) x$meetingtime)
  print(summary(meetingtimes))
  meantau <- mean(meetingtimes)
  mediantau <- median(meetingtimes)
  k <- floor(as.numeric(quantile(meetingtimes, probs = 0.90)) + 1)
  K <- floor(10*k)
  cat("found k =", k, "and K =", K, "\n")
  clusterExport(cl, "K")
  c_chains_ <- parLapply(cl, 1:nsamples, function(x) f(x, lambda, K))
  mean_estimators <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
    H_bar(c_chains_[[irep]], h = function(x) x[1:p], k = k, K = K)
  }
  sum_estimators <- colSums(mean_estimators)
  sum_squareestimators <- colSums(mean_estimators^2)
  df <- rbind(df, data.frame(ilambda = rep(ilambda, p), lambda = rep(lambda, p), component = 1:p, k = rep(k, p),
                             K = rep(K, p), nsamples = rep(nsamples, p),
                             meantau = meantau,
                             mediantau = mediantau,
                             sum_est = matrix(sum_estimators, nrow = p),
                             sum_squareest = matrix(sum_squareestimators, nrow = p)))
  save(df, nsamples, lambdas, file = "diabetes.cchains.p10.RData")
}
load(file = "diabetes.cchains.p10.RData")

stopCluster(cl)

head(df)
