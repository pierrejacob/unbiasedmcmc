# load packages
library(debiasedmcmc)
# library(lars)
setmytheme()
rm(list = ls())
set.seed(21)
data(diabetes);

# save(diabetes, file = "~/Dropbox/PolyaGamma/code/debiasedmcmc/data/diabetes.RData")
X <- scale(diabetes$x2);
Y <- matrix(scale(diabetes$y), ncol = 1);
f <- function(x, lambda){
  pb <- get_blasso(Y, X, lambda);
  return(coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit))
}


library(parallel)
RNGkind("L'Ecuyer-CMRG")
mc.cores <- 6

nsamples <- 50

f(1,1)
results <- mclapply(1:nsamples,
                    function(x) rnorm(100),
                    mc.cores=mc.cores)

# sapply(results, function(x) x$meetingtime) %>% summary
###
library(parallel)
cl <- makeCluster(6, type = "FORK")
clusterSetRNGStream(cl, 12)

f(1, 1)
nsamples <- 50
c_chains_ <- parLapply(cl, 1:nsamples, function(x) f(x, 1))
stopCluster(cl)



clusterSetRNGStream(cl, 12)


cl <- makeCluster(detectCores()-2)
clusterSetRNGStream(cl, 12)

clusterEvalQ(cl, library("debiasedmcmc"))
clusterEvalQ(cl, library("lars"))
clusterEvalQ(cl, {data(diabetes);
  X <- scale(diabetes$x2);
  Y <- matrix(scale(diabetes$y), ncol = 1);
})

clusterCall(cl, runif, 3)
parLapply(cl, 1:15, fun = function(x) runif(3))

## now try coupled chains

f <- function(x, lambda){
  pb <- get_blasso(Y, X, lambda);
  return(coupled_chains(pb$gibbs_kernel, pb$coupled_gibbs_kernel, pb$rinit))
}
data(diabetes);
X <- scale(diabetes$x2);
Y <- matrix(scale(diabetes$y), ncol = 1);

clusterExport(cl, "f")
nsamples <- 50
proc.time()
for (i in 1:nsamples){
  e <- f(x, 1)
}
proc.time()
c_chains_ <- parLapply(cl, 1:nsamples, function(x) f(x, 1))
proc.time()
stopCluster(cl)

meetingtimes <- sapply(c_chains_, function(x) x$meetingtime)
registerDoParallel(cores = detectCores())
mean_estimators <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
  H_bar(c_chains_[[irep]], h = function(x) x[1:5], k = 0, K = 1)
}
colMeans(mean_estimators)
print(summary(meetingtimes))
