# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
# Plummer example
# nhpv considered as Y
nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
y1 <- matrix(nhpv, nrow = 1)
# Npart is put in the parameters
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
           143, 229, 696, 93)
J <- 13
# posterior is beta in each study, with parameters
posterior_phi_alpha <- 1 + nhpv
posterior_phi_beta <- 1 + Npart - nhpv

sample_module1 <- function(nsamples){
  theta1s <- matrix(nrow = nsamples, ncol = J)
  for (j in 1:J){
    theta1s[,j] <- rbeta(nsamples, shape1 = posterior_phi_alpha[j], shape2 = posterior_phi_beta[j])
  }
  return(theta1s)
}

###
# For module 2, ncases considered data
ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
y2 <- matrix(ncases, nrow = 1)
# Npop considered parameters
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
          26751, 75815, 150302, 354993, 3683043, 507218)
Npop_normalized <- log(10**(-3) * Npop)
# Find parameters given Y2
hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = sqrt(1000))
dprior2 <- function(theta2, hyper2){
  return(sum(dnorm(theta2, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)))
}

cppFunction("double plummer_module2_loglikelihood_(NumericVector theta1, NumericVector theta2,
                                           NumericVector ncases, NumericVector Npop_normalized){
  double eval = 0;
  double mu, logmu;
  for (int j = 0; j < 13; j++){
    logmu = theta2(0) + theta1(j) * theta2(1) + Npop_normalized(j);
    mu = exp(logmu);
    eval += ncases(j) * logmu - mu;
  }
  return eval;
}
")

get_kernels <- function(theta1, Sigma_proposal){
  target <- function(x) plummer_module2_loglikelihood_(theta1, x, ncases, Npop_normalized) + dprior2(x, hyper2)
  ##
  # Markov kernel of the chain
  single_kernel <- function(chain_state){
    proposal_value <- chain_state + fast_rmvnorm(1, rep(0, dimension), Sigma_proposal)[1,]
    proposal_pdf <- target(proposal_value)
    current_pdf <- target(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(proposal_value)
    } else {
      return(chain_state)
    }
  }

  # Markov kernel of the coupled chain
  # Sigma <- diag(sd_proposal^2, dimension, dimension)
  Sigma_chol <- t(chol(Sigma_proposal))
  Sigma_chol_inv <- t(solve(chol(Sigma_proposal)))
  coupled_kernel <- function(chain_state1, chain_state2){
    distance_ <- mean((chain_state1 - chain_state2)^2)
    if (distance_ > 5){
      proposal_value <- gaussian_opt_transport(1, chain_state1, chain_state2, Sigma_chol, Sigma_chol, Sigma_chol_inv, Sigma_chol_inv)[[1]]
      proposal1 <- proposal_value[,1]
      proposal2 <- proposal_value[,2]
    } else {
      proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
      proposal1 <- proposal_value[,1]
      proposal2 <- proposal_value[,2]
    }
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- target(proposal2)
    current_pdf1 <- target(chain_state1)
    current_pdf2 <- target(chain_state2)
    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
    }
    if (accept2){
      chain_state2 <- proposal2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  return(list(target = target, coupled_kernel = coupled_kernel, single_kernel = single_kernel))
}

# theta1s <- sample_module1(1)
# theta1 <- theta1s[1,]
# dimension <- 2
# sd_proposal <- 1
# kernels <- get_kernels(theta1, sd_proposal)
# single_kernel <- kernels$single_kernel
# coupled_kernel <- kernels$coupled_kernel
# target <- kernels$target
#
# rinit <- function() fast_rmvnorm(1, mean = rep(0, dimension, dimension), diag(1, dimension, dimension))[1,]
# ### test of the coupled kernel
# niterations <- 10000
# current_value1 <- rinit()
# current_value2 <- rinit()
# current_pdf1 <- target(current_value1)
# current_pdf2 <- target(current_value2)
# chain1 <- matrix(ncol=dimension, nrow=niterations)
# chain2 <- matrix(ncol=dimension, nrow=niterations)
# chain1[1,] <- current_value1
# chain2[1,] <- current_value2
# for (t in 2:niterations){
#   current_value <- coupled_kernel(current_value1, current_value2)
#   current_value1 <- current_value$chain_state1
#   current_value2 <- current_value$chain_state2
#   chain1[t,] <- current_value1
#   chain2[t,] <- current_value2
# }
#
# distances_ <- sapply(1:niterations, function(t) mean((chain1[t,] - chain2[t,])^2))
# plot(1:niterations, distances_, type = "l")
#
# it <- floor(seq(from = 1, to  = niterations, length.out = 1000))
# qplot(x=it, y=chain1[it,1], geom = "line") + ylab("X") + xlab("iteration") + geom_line(aes(y = chain2[it,1]), colour = "red")
# qplot(x=it, y=chain1[it,2], geom = "line") + ylab("X") + xlab("iteration") + geom_line(aes(y = chain2[it,2]), colour = "red")
#
# hist(chain1[100:niterations,1])
# hist(chain1[100:niterations,2])

##
dimension <- 2
Sigma_proposal <- diag(c(1,4), dimension, dimension)
nsamples <- 500
theta1s <- sample_module1(nsamples)
rinit <- function() fast_rmvnorm(1, mean = c(0, 0), diag(1, dimension, dimension))[1,]

c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  theta1 <- theta1s[irep,]
  kernels <- get_kernels(theta1, Sigma_proposal)
  single_kernel <- kernels$single_kernel
  coupled_kernel <- kernels$coupled_kernel
  coupled_chains(single_kernel, coupled_kernel, rinit)
}

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)

# ##
K <- 1000
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  theta1 <- theta1s[irep,]
  kernels <- get_kernels(theta1, Sigma_proposal)
  single_kernel <- kernels$single_kernel
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}
#
k <- 500

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], k = k, K = K)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], h = function(x) x^2, k = k, K = K)
}

est_mean <- rep(0, dimension)
est_var <- rep(0, dimension)
for (component in 1:dimension){
  estimators <- sapply(mean_estimators, function(x) x[component])
  est_mean[component] <- mean(estimators)
  cat("estimated mean: ", est_mean[component], "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
  s_estimators <- sapply(square_estimators, function(x) x[component])
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
  cat("estimated variance: ", est_var[component], "\n")
}

Sigma_proposal <- diag(est_var, dimension, dimension)
# Sigma_proposal[1,2] <- Sigma_proposal[2,1] <- -0.3
nsamples <- 5000
theta1s <- sample_module1(nsamples)
dimension <- 2
rinit <- function() fast_rmvnorm(1, mean = est_mean, Sigma_proposal)[1,]

c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  theta1 <- theta1s[irep,]
  kernels <- get_kernels(theta1, Sigma_proposal)
  single_kernel <- kernels$single_kernel
  coupled_kernel <- kernels$coupled_kernel
  coupled_chains(single_kernel, coupled_kernel, rinit)
}

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)
k <- 2*floor(as.numeric(quantile(meetingtime, probs = 0.95)))
K <- 10*k
# # ind_ <- which(sapply(c_chains_, function(x) x$meetingtime > 1000))
# chain1 <- c_chains_[[1]]$samples1
# chain2 <- c_chains_[[1]]$samples2
#
# qplot(x=0:(nrow(chain1)-1), y=chain1[,1], geom = "line") + ylab("X") + xlab("iteration") +
#   geom_line(aes(x = 1:nrow(chain2), y = chain2[,1]), colour = "red")
#
# qplot(x=0:(nrow(chain1)-1), y=chain1[,2], geom = "line") + ylab("X") + xlab("iteration") +
#   geom_line(aes(x = 1:nrow(chain2), y = chain2[,2]), colour = "red")


# ##
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  theta1 <- theta1s[irep,]
  kernels <- get_kernels(theta1, Sigma_proposal)
  single_kernel <- kernels$single_kernel
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}

mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], k = k, K = K)
}

square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], h = function(x) x^2, k = k, K = K)
}

cross_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], h = function(x) x[1] * x[2], k = k, K = K)
}

est_mean <- rep(0, dimension)
est_var <- rep(0, dimension)

for (component in 1:dimension){
  estimators <- sapply(mean_estimators, function(x) x[component])
  cat("estimated mean: ", mean(estimators), "+/- ", 2*sd(estimators)/sqrt(nsamples), "\n")
  s_estimators <- sapply(square_estimators, function(x) x[component])
  cat("estimated second moment: ", mean(s_estimators), "+/- ", 2*sd(s_estimators)/sqrt(nsamples), "\n")
  cat("estimated variance: ", mean(s_estimators) - mean(estimators)^2, "\n")
  est_mean[component] <- mean(estimators)
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
}
c_estimators <- sapply(cross_estimators, function(x) x[1])
cat("estimated covariance: ", mean(c_estimators) - prod(est_mean), "\n")

# histogram
histogram1 <- histogram_c_chains(c_chains_continued_, 1, k, K, nclass = 100)
plot_histogram(histogram1, with_bar = F) + xlab(expression(theta[1])) + ylab("density")
# histogram <- histogram_c_chains(c_chains_continued_, 1, k, K, breaks = seq(from = -1.8, to = -1.3, length.out = 100))
# g <- plot_histogram(histogram)
# g

histogram2 <- histogram_c_chains(c_chains_continued_, 2, k, K, nclass = 100)
plot_histogram(histogram2, with_bar = F) + xlab(expression(theta[2])) + ylab("density")
# histogram <- histogram_c_chains(c_chains_continued_, 2, k, K, breaks = seq(from = 10, to = 15, length.out = 100))
# g <- plot_histogram(histogram)
# g
