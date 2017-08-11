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
    proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
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

nsamples <- 1000
theta1s <- sample_module1(nsamples)
theta1hat <- colMeans(theta1s)
##
dimension <- 2
Sigma_proposal <- diag(c(1,1), dimension, dimension)
rinit <- function() fast_rmvnorm(1, mean = c(0, 0), diag(1, dimension, dimension))[1,]

filename <- "plummer.tuning.RData"
kernels <- get_kernels(theta1hat, Sigma_proposal)
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  # theta1 <- theta1s[irep,]
  single_kernel <- kernels$single_kernel
  coupled_kernel <- kernels$coupled_kernel
  coupled_chains(single_kernel, coupled_kernel, rinit)
}
save(c_chains_, file = filename)
load(file = filename)
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
# hist(meetingtime)

k <- as.numeric(floor(quantile(meetingtime, probs = 0.95)))
K <- 2*k

#
# K <- 1000
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  single_kernel <- kernels$single_kernel
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}
save(nsamples, k, K, c_chains_, c_chains_continued_, file = filename)
load(file = filename)
#
# k <- 500

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
print(Sigma_proposal)


filename <- "plummer.results.RData"
nsamples <- 10000

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
save(nsamples, c_chains_, file = filename)
load(filename)
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
# sum(sapply(c_chains_2, function(x) x$iteration))
# hist(meetingtime)

# x <- as.numeric(names(table(meetingtime)))
# y <- as.numeric(table(meetingtime)) / nsamples
# g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
# g
g <- qplot(x = meetingtime, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab("meeting time") + ylab("proportion")
g
ggsave(filename = "plummer.meetingtimes.pdf", plot = g, width = 5, height = 5)


k <- as.numeric(floor(quantile(meetingtime, probs = 0.95)))
K <- 10*k
# nsamples <- 1000
# ##
c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  theta1 <- theta1s[irep,]
  kernels <- get_kernels(theta1, Sigma_proposal)
  single_kernel <- kernels$single_kernel
  continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
}
save(k, K, nsamples, c_chains_, c_chains_continued_, file = filename)

# save(c_chains_, c_chains_continued_, c_chains_2, c_chains_continued_2, file = filename)
load(filename)
sum(sapply(c_chains_continued_, function(x) x$iteration))


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

# # histogram
# histogram1 <- histogram_c_chains(c_chains_continued_2, 1, k, K, nclass = 100)
# g1 <- plot_histogram(histogram1, with_bar = F) + xlab(expression(theta[2.1])) + ylab("density")
# g1
# ggsave("plummer.histogram1.pdf", plot = g1, width = 7, height = 7)
#
# # g
#
# histogram2 <- histogram_c_chains(c_chains_continued_2, 2, k, K, nclass = 100)
# g2 <- plot_histogram(histogram2, with_bar = F) + xlab(expression(theta[2.2])) + ylab("density")
# g2
# ggsave("plummer.histogram2.pdf", plot = g2, width = 7, height = 7)
# # g

### exact cut distribution from tedious parallel MCMC
niterations <- 1000
theta2s <- foreach(itheta = 1:nrow(theta1s), .combine = rbind) %dorng% {
  theta1 <- theta1s[itheta,]
  kernels <-  get_kernels(theta1, Sigma_proposal)
  chain <- rinit()
  for (iter in 1:niterations){
    chain <- kernels$single_kernel(chain)
  }
  chain
}
save(theta2s, file = "plummer.mcmc.RData")
# load(file = "plummer.mcmc.RData")

# hist_mcmc <- hist(theta2s[,1], breaks = histogram1$breaks, plot = F)
# g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[2.1])) + ylab("density")
# g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
# g1
# ggsave("plummer.histogram1.pdf", plot = g1, width = 7, height = 7)
#
# hist_mcmc <- hist(theta2s[,2], breaks = histogram2$breaks, plot = F)
# g2 <- plot_histogram(histogram2, with_bar = T) + xlab(expression(theta[2.2])) + ylab("density")
# g2 <- g2 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
# g2
# ggsave("plummer.histogram2.pdf", plot = g2, width = 7, height = 7)
