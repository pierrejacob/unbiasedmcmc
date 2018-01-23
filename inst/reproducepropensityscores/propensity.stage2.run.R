# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(8)
registerDoParallel(cores = detectCores())
#
load("propensity.data.RData")
Y <- logistic_setting$Y
X <- logistic_setting$X
p <- ncol(X)
n <- nrow(X)

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

## Stage 2 given plug-in estimate of first stage
# first, come up with plug-in estimate
filename <- "propensity.stage1.c_chains.RData"
load(file = filename)
nsamples <- length(c_chains_continued_)
mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], k = k, K = K)
}
theta1_hat <- rep(0, p)
for (component in 1:p){
  theta1_hat[component] <- mean(sapply(mean_estimators, function(x) x[component]))
}
theta1_hat

# k <- 20
# K <- 5*k
# mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
#   H_bar(c_chains_continued_[[irep]], k = k, K = K)
# }
# square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
#   H_bar(c_chains_continued_[[irep]], h = function(x) x^2, k = k, K = K)
# }
# est_mean <- rep(0, p)
# est_var <- rep(0, p)
# for (component in 1:p){
#   cat("component ", component, "\n")
#   estimators <- sapply(mean_estimators, function(x) x[component])
#   cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
#   s_estimators <- sapply(square_estimators, function(x) x[component])
#   cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
#   est_mean[component] <- mean(estimators)
#   est_var[component] <- mean(s_estimators) - est_mean[component]^2
# }

# beta_hat <- est_mean

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

logitlink <- apply(X = X, MARGIN = 1, FUN = function(row)
  0.6 * row[1] + 0.5 * exp(row[2] - 1) + 0.4 * row[3] +
    0.3 * exp(row[4] - 1) + 0.2 * abs(row[5]) + 0.1 * abs(row[6]))
second_Y <- rbinom(n, 1, expit(logitlink))
filename <- "propensity.second_data.RData"
save(second_Y, file = filename)

#
second_b <- matrix(0, nrow = 6, ncol = 1)
second_B <- diag(50, 6, 6)
second_B[1,1] <- 800
second_p <- 6
#
# #
# second_X <- get_second_X(beta_hat)
# second_logistic_setting <- logistic_precomputation(second_Y, second_X, second_b, second_B)

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
  coupled_kernel <- function(chain_state1, chain_state2){
    ws <- sample_w(chain_state1, chain_state2, second_logistic_setting$X)
    betas <- sample_beta(ws$w1, ws$w2, second_logistic_setting)
    if (all(ws$w1 == ws$w2)){
      betas$beta2 <- betas$beta1
    }
    return(list(chain_state1=cbind(betas$beta1),
                chain_state2=cbind(betas$beta2)))
  }
  rinit <- function() fast_rmvnorm(1, rep(0, second_p), diag(1, second_p, second_p))
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit))
}

# pb <- get_setting(theta1_hat)
# niterations <- 10000
# betas <- matrix(0, ncol=second_p, nrow=niterations)
# beta <- pb$rinit()
# for (i in 1:niterations){
#   beta <- pb$single_kernel(beta)
#   betas[i,] <- beta
# }
#
# matplot(betas, type = "l")
# hist(betas[,2], nclass = 100)

## now the cut approach
sapply(c_chains_continued_, function(x) x$meeting) %>% summary
sapply(c_chains_continued_, function(x) x$iteration) %>% summary

# coupled_chains(pb$single_kernel, pb$coupled_kernel, pb$rinit)

nsamples <- length(c_chains_continued_)
second_c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  c1_ <- c_chains_continued_[[irep]]
  beta1 <- c1_$samples1[k,]
  second_stage <- get_setting(beta1)
  coupled_chains(second_stage$single_kernel, second_stage$coupled_kernel, second_stage$rinit)
}

sapply(second_c_chains_, function(x) x$meeting) %>% summary

k2 <- as.numeric(floor(quantile(sapply(second_c_chains_, function(x) x$meeting), 0.95)))
K2 <- 10*k2
second_c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  c1_ <- c_chains_continued_[[irep]]
  beta1 <- c1_$samples1[k,]
  second_stage <- get_setting(beta1)
  coupled_chains(second_stage$single_kernel, second_stage$coupled_kernel, second_stage$rinit, K = K2)
}
save(nsamples, k2, K2, second_c_chains_, file = "propensity.stage2.c_chains.RData")

sapply(second_c_chains_, function(x) x$meeting) %>% summary
sapply(second_c_chains_, function(x) x$iteration) %>% summary

meetingtime <- sapply(second_c_chains_, function(x) x$meetingtime)
x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g

histogram1 <- histogram_c_chains(second_c_chains_, 2, k2, K2, nclass = 20)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta)) + ylab("density")
g1
# compared to plug-in
# g1 + geom_density(data=data.frame(x = betas[,2]), aes(x = x, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, y = ..density..))
