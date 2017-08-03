
# Setup -------------------------------------------------------------------
rm(list = ls())

library(debiasedmcmc)
library(corpcor)

setmytheme()
set.seed(22)
registerDoParallel(cores = detectCores())


# Draw input data ---------------------------------------------------------
n <- 1000
p <- 40

# covariates
meanX <- rep(0, p)
sigmaX <- diag(1, nrow = p, ncol = p)
X <- fast_rmvnorm(n, meanX, sigmaX)

# outcome
beta_star <- 1:p / p
logitprobs <- apply(X = X, MARGIN = 1, FUN = function(row) sum(row * beta_star))
Y <- rbinom(n, 1, expit(logitprobs))

# prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
logistic_setting <- logistic_precomputation(Y, X, b, B)


# Run PGG and build rinit functions --------------------------------------

# PGG estimation
gibbs_iters = 1000
betas_pg <- pg_gibbs(gibbs_iters, logistic_setting)
matplot(betas_pg, type = "l", lty = 1)

# build rinits
m_hat <- colMeans(betas_pg[100:gibbs_iters,])
cov_hat <- cov.shrink(betas_pg[100:gibbs_iters,], verbose=FALSE)

rinit_zeroes <- function() matrix(rep(0,p),p,1)
rinit_bad_prior <- function() t(fast_rmvnorm(1, mean = matrix(10, nrow=p, ncol=1), covariance = B))
rinit_prior <- function() t(fast_rmvnorm(1, mean = b, covariance = B))
rinit_post <- function() t(fast_rmvnorm(1, mean = m_hat, covariance = cov_hat))

# Markov kernel of a single chain
single_kernel <- function(chain_state, logistic_setting){
  zs <- abs(xbeta(logistic_setting$X, t(chain_state)))
  w <- rpg(logistic_setting$n, h=1, z=zs)
  res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(chain_state)
}

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2, logistic_setting){
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  return(list(chain_state1=cbind(betas$beta1),
              chain_state2=cbind(betas$beta2)))
}

# Define scalable coupled kernel -------------------------

make_cond_normal_data <- function(res1,res2){
  mu1 <- res1$m
  mu2 <- res2$m
  Sigma1 = tcrossprod(res1$Cholesky)
  Sigma2 = tcrossprod(res2$Cholesky)

  cond_normal_data = vector(p,mode='list')
  for(j in 1:p){
    in_idx <- j
    out_idx <- setdiff(1:p,j)
    m1_proj <- Sigma1[in_idx,out_idx] %*% solve(Sigma1[out_idx,out_idx])
    m2_proj <- Sigma2[in_idx,out_idx] %*% solve(Sigma2[out_idx,out_idx])
    S1_cond <- Sigma1[in_idx,in_idx] - m1_proj %*% Sigma1[out_idx,in_idx]
    S2_cond <- Sigma2[in_idx,in_idx] - m2_proj %*% Sigma2[out_idx,in_idx]
    Cholesky1_cond <- sqrt(S1_cond)
    Cholesky2_cond <- sqrt(S2_cond)
    Cholesky1_inv_cond <- 1/Cholesky1_cond
    Cholesky2_inv_cond <- 1/Cholesky2_cond
    cond_normal_data[[j]] = list(in_idx=in_idx,out_idx=out_idx,
                                 mu1=mu1, mu2=mu2, m1_proj=m1_proj, m2_proj=m2_proj,
                                 Cholesky1=Cholesky1_cond, Cholesky2=Cholesky2_cond,
                                 Cholesky1_inv=Cholesky1_inv_cond, Cholesky2_inv=Cholesky2_inv_cond)
  }
  return(cond_normal_data)
}

max_coupling_cycle <- function(beta1,beta2,cond_normal_data){
  for(j in 1:p){
    d <- cond_normal_data[[j]]
    out_idx <- d$out_idx
    mu1 <- d$mu1
    mu2 <- d$mu2
    m1 <- mu1[d$in_idx] + d$m1_proj %*% (beta1[out_idx]-mu1[out_idx])
    m2 <- mu2[d$in_idx] + d$m2_proj %*% (beta2[out_idx]-mu2[out_idx])
    x <- debiasedmcmc:::gaussian_max_coupling_cholesky(m1, m2, d$Cholesky1, d$Cholesky2, d$Cholesky1_inv, d$Cholesky2_inv)
    beta1[j] <- x[,1]
    beta2[j] <- x[,2]
  }
  return(list(beta1=beta1, beta2=beta2))
}

sample_beta_scalable <- function(w1, w2, beta1, beta2, logistic_setting, mode="scalable"){
  X <- logistic_setting$X
  KTkappaplusinvBtimesb <- logistic_setting$KTkappaplusinvBtimesb
  invB <- logistic_setting$invB
  res1 <- m_and_sigma(w1, X, invB, KTkappaplusinvBtimesb)
  res2 <- m_and_sigma(w2, X, invB, KTkappaplusinvBtimesb)

  if(mode=='scalable'){
    p <- logistic_setting$p
    #nd <- floor(1+log(p))
    nd <- 1
    cond_normal_data <- make_cond_normal_data(res1,res2)
    for(k in 1:nd){
      x <- max_coupling_cycle(beta1,beta2,cond_normal_data)
      beta1 <- x$beta1
      beta2 <- x$beta2
    }
  } else {
    stop('invalid coupling method')
  }
  return(list(beta1=beta1, beta2=beta2))
}

# Scalable Markov kernel of the coupled chain
scalable_coupled_kernel <- function(chain_state1, chain_state2, logistic_setting){
 ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
 betas <- sample_beta_scalable(ws$w1, ws$w2, chain_state1, chain_state2, logistic_setting)
 if (all(ws$w1 == ws$w2)){
   betas$beta2 <- betas$beta1
 }
 return(list(chain_state1=cbind(betas$beta1),
             chain_state2=cbind(betas$beta2)))
}




# Test of the coupled kernel ----------------------------------------------
niterations <- 1000
current_value1 <- rinit_post()
current_value2 <- rinit_post()

chain1 <- matrix(ncol=p, nrow=niterations)
chain2 <- matrix(ncol=p, nrow=niterations)
chain1[1,] <- current_value1
chain2[1,] <- current_value2
for (t in 2:niterations){
  current_value <- coupled_kernel(current_value1, current_value2, logistic_setting)
  chain1[t,] <- current_value1 <- current_value$chain_state1
  chain2[t,] <- current_value2 <- current_value$chain_state2
}

distances_ <- sapply(1:niterations, function(t) mean((chain1[t,] - chain2[t,])^2))
plot(1:niterations, distances_, type = "l")


# Test of the scalable kernel ----------------------------------------------
niterations <- 1000
current_value1 <- rinit_post()
current_value2 <- rinit_post()

chain1 <- matrix(ncol=p, nrow=niterations)
chain2 <- matrix(ncol=p, nrow=niterations)
chain1[1,] <- current_value1
chain2[1,] <- current_value2
for (t in 2:niterations){
  current_value <- scalable_coupled_kernel(current_value1, current_value2, logistic_setting)
  chain1[t,] <- current_value1 <- current_value$chain_state1
  chain2[t,] <- current_value2 <- current_value$chain_state2
}

distances_ <- sapply(1:niterations, function(t) mean((chain1[t,] - chain2[t,])^2))
plot(1:niterations, distances_, type = "l")







# Find optimal k and K ----------------------------------------------------

R <- 1000
c_chains_ <-  foreach(irep = 1:R) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit_post, logistic_setting)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)

# Select optimal k
asympt_ineffic_fn <- function(m,observed_tau){
  e_max <- mean(pmax(observed_tau,m))
  asymptotic_ineffic <- e_max * (1 + 2*e_max - m)
  return(asymptotic_ineffic)
}
# would be good to estimate these variances directly as well, as in the latter part of v1 of the paper

observed_tau <- meetingtime
k_vec <- 1:max(observed_tau)
ae_vec <- sapply(k_vec,asympt_ineffic_fn, observed_tau=observed_tau)
df <- data.frame(k=k_vec,ae=ae_vec)
k_opt <- k_vec[which.min(ae_vec)]
df <- df[df$k < k_opt*4,]
ggplot(df,aes(x=k,y=ae)) +
  geom_line() + geom_vline(xintercept=k_opt,col='red') +
  ylab('Asymptotic Inefficiency')

K_opt <- k_opt * 10

k_opt
K_opt

# Compare unbiased estimator with truth ---------------------------------

c_chains_ <- foreach(chain = c_chains_) %dorng% {
  continue_coupled_chains(chain, single_kernel, K_opt, logistic_setting)
}

mean_estimators <-  foreach(chain = c_chains_) %dorng% {
  H_bar(chain, k=k_opt, K=K_opt)
}
square_estimators <-  foreach(chain = c_chains_) %dorng% {
  H_bar(chain, h = function(x) x^2, k=k_opt, K=K_opt)
}

est_mean <- rep(0, p)
est_var <- rep(0, p)
for (component in 1:p){
  estimators <- sapply(mean_estimators, function(x) x[component])
  s_estimators <- sapply(square_estimators, function(x) x[component])
  est_mean[component] <- mean(estimators)
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
  cat("\ncomponent ", component, "\n")
  cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(length(estimators)), "\n")
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(length(estimators)), "\n")
  cat("estimated variance: ", mean(s_estimators) - mean(estimators)^2, "\n")
}



compar.df <- data.frame(method="unbiased",component=paste0("estimate",1:p),
                        truth=beta_star,est_mean=est_mean,est_var=est_var)
compar.df.mcmc <- data.frame(method="mcmc",component=paste0("estimate",1:p),
                             truth=beta_star,est_mean=m_hat,est_var=est_var)
compar.df <- rbind(compar.df,compar.df.mcmc)
compar.df <- compar.df %>% mutate(ymin=est_mean-2*sqrt(est_var), ymax=est_mean+2*sqrt(est_var))

g <- ggplot(compar.df, aes(x = truth, y = truth)) + geom_line() + geom_point() + theme_bw()
g <- g + geom_errorbar(aes(ymin = ymin, ymax = ymax, group = method, colour = method)) +
  scale_color_manual(values = c('black', "orange"))
g




