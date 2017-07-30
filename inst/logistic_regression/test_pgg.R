
# Setup -------------------------------------------------------------------
library(debiasedmcmc)
library(corpcor)

rm(list = ls())

setmytheme()
set.seed(21)
registerDoParallel(cores = detectCores())


# Draw input data ---------------------------------------------------------
n <- 1000
p <- 25

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
  chain1[t,] <- current_value$chain_state1
  chain2[t,] <- current_value$chain_state2
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

# Compare unbiased estimator with truth ---------------------------------

c_chains_

x <- c_chains_ <- foreach(irep = 1:R) %dorng% {
}


x <- continue_coupled_chains(c_chains_[[1]], single_kernel, K=K_opt)



estimates.long <- estimates %>% gather(component, value, estimate.1:estimate.5)
# summary <- estimates.long %>% group_by(component) %>% summarize(iter = mean(iteration), m = mean(value), sd = sd(value)/sqrt(nrep))
# summary$method <- "unbiased"
#
# estimates.m <- estimates_m %>% gather(component, value, estimate.1:estimate.5)
# summary_m <- estimates.m %>% group_by(component) %>% summarize(iter = mean(iteration), m = mean(value), sd = sd(value)/sqrt(nrep))
# summary_m$method <- "unbiased+m"
#
# mcmc <- data.frame(truth = colMeans(betas_pg), component = paste0("estimate.", 1:p))
# compar.df <- merge(mcmc, summary, by = "component")
# compar.df_m <- merge(mcmc, summary_m, by="component")
# compar.df <- rbind(compar.df, compar.df_m)
#
# compar.df <- compar.df %>% mutate(ymin = m - 2*sd, ymax = m + 2*sd)
#
# g <- ggplot(compar.df, aes(x = truth, y = truth)) + geom_line() + geom_point() + theme_bw()
# g <- g + geom_errorbar(aes(ymin = ymin, ymax = ymax, group = method, colour = method)) + scale_color_manual(values = c('black', "orange"))
# g






