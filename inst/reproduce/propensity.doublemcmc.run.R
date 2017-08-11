# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores())

cppFunction("
NumericVector beta2e_(const NumericVector & beta, const NumericMatrix & C){
            double betatimescovar;
            int nbeta = beta.size();
            int nobservations = C.nrow();
            NumericVector e(nobservations);
            for (int i = 0; i < nobservations; i++){
            betatimescovar = beta(0);
            for (int j = 0; j < nbeta-1; j++){
            betatimescovar += beta(j+1) * C(i,j);
            }
            e(i) = 1. / (1. + exp(-betatimescovar));
            }
            return e;
            }")

# function that takes a vector of doubles
# compute the quintiles
# and return a vector of indicators of quintile memberships,
# ie with values in {1,2,3,4,5}.
cppFunction("
IntegerVector cut_in_fifth_(const NumericVector & x){
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  NumericVector q = NumericVector::create(0.2, 0.4, 0.6, 0.8);
  int s = x.size();
  int nqs = q.size();
  NumericVector quintiles(nqs);
  for (int i = 0; i < nqs; i ++){
    quintiles(i) = y[(int) (s * q[i])];
  }
  IntegerVector indicators(s);
  int indicator;
  for (int i = 0; i < s; i ++){
    indicator = 0;
    while ((x(i) >= quintiles(indicator)) && (indicator < (nqs - 1))){
      indicator ++;
    }
    if (x(i) >= quintiles(indicator)){ indicator ++; }
    indicators(i) = indicator + 1; // so that the result is between 1 and 5.
  }
  return indicators;
}")

#

load(file = "propensity.data.RData")
load(file = "propensity.second_data.RData")
p <- logistic_setting$p
X <- logistic_setting$X
Y <- logistic_setting$Y
n <- nrow(X)
#
second_b <- matrix(0, nrow = 6, ncol = 1)
second_B <- diag(50, 6, 6)
second_B[1,1] <- 800
second_p <- 6
#

single_kernel <- function(chain_state){
  beta <- chain_state[1:p]
  zs <- abs(xbeta(X, beta))
  w <- rpg(n, h=1, z=zs)
  res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  beta <- fast_rmvnorm_chol(1, res$m, res$Cholesky)
  return(beta)
}

get_second_X <- function(beta){
  # propensity scores
  e <- beta2e(beta, X)
  indicators <- cut_in_fifth(e)
  # now regression of Z on intercept, Y, and indicators
  second_X <- cbind(rep(1, n), Y, matrix(0, nrow = n, ncol = 4))
  for (i in 2:5){
    second_X[,1+i] <- (indicators == i)
  }
  return(second_X)
}

get_setting <- function(theta1){
  second_X <- get_second_X(theta1)
  second_logistic_setting <- logistic_precomputation(second_Y, second_X, second_b, second_B)
  single_kernel <- function(chain_state){
    beta <- chain_state[1:second_p]
    zs <- abs(xbeta(second_X, beta))
    w <- rpg(n, h=1, z=zs)
    res <- m_and_sigma(w, second_X, second_logistic_setting$invB, second_logistic_setting$KTkappaplusinvBtimesb)
    beta <- fast_rmvnorm_chol(1, res$m, res$Cholesky)
    return(beta)
  }
  return(single_kernel)
}

niterations1 <- 1e5
# betas <- matrix(0, ncol=p, nrow=niterations1)
# beta <- rep(0, p)
# for (i in 1:niterations1){
#   beta <- single_kernel(beta)
#   betas[i,] <- beta
# }
# save(niterations1, betas, file = "propensity.doublemcmc.RData")
load("propensity.doublemcmc.RData")
matplot(betas[1:100,], type = "l")

# beta1 <- rep(0, p)
# for (i in 1:niterations){
#   beta1 <- single_kernel(beta1)
# }
# beta1 <- betas[niterations,]
# kernel <- get_setting(beta1)
# betas2 <- matrix(0, ncol=second_p, nrow=niterations)
# beta2 <- rep(0, second_p)
# for (i in 1:niterations){
#   beta2 <- kernel(beta2)
#   betas2[i,] <- beta2
# }
# matplot(betas2, type = "l")

nsamples <- 10000
# get sub sample post burnin of 1e4
subsample <- floor(seq(from = 1e4, to = niterations1, length.out = nsamples))
betas1 <- betas[subsample,,drop=FALSE]

niterations <- 100
double_mcmc <-  foreach(irep = 1:nsamples, .combine = rbind) %dorng% {
  beta1 <- betas1[irep,]
  kernel <- get_setting(beta1)
  beta2 <- rep(0, second_p)
  for (i in 1:niterations){
    beta2 <- kernel(beta2)
  }
  beta2
}
save(niterations1, betas, subsample, betas1, nsamples, niterations, double_mcmc, file = "propensity.doublemcmc.RData")
hist(double_mcmc[,2])
