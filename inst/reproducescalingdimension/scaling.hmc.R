library(unbiasedmcmc)
library(doParallel)
library(doRNG)
# setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores()-2)

# leapfrog specifically for MVN target with covariance matrix 'variance', precision matrix 'precision', zero mean
cppFunction('
            Rcpp::List leapfrog_c(Eigen::VectorXd & x, Eigen::VectorXd & v, const Eigen::MatrixXd & precision,
            const Eigen::MatrixXd & variance, int nsteps, double stepsize){
            v = v + ((-precision * x ).array() * stepsize / 2.).matrix();
            for (int step = 0; step < nsteps; step ++){
            x = x + ((variance * v).array() * stepsize).matrix();
            if (step != nsteps){
            v = v + ((-precision * x).array() * stepsize).matrix();
            }
            }
            v = v + ((-precision * x).array() * stepsize / 2.).matrix();
            return(Rcpp::List::create(Rcpp::Named("x") = x, Rcpp::Named("v") = v));
            }
            ', depends = "RcppEigen")

cppFunction('
            Eigen::VectorXd gradlogtarget_c(const Eigen::VectorXd & x, const Eigen::VectorXd & mean, const Eigen::MatrixXd & precision){
            return(- precision * (x - mean));
            }
            ', depends = "RcppEigen")


# One step of HMC
hmc_single_kernel <- function(chain_state, current_pdf, nsteps, stepsize, pb){
  dimension_ <- length(chain_state)
  current_v <- fast_rmvnorm(1, rep(0, dimension_), pb$precision_pi)
  leapfrog_result <- leapfrog_c(chain_state, current_v, pb$precision_pi, pb$Sigma_pi, nsteps, stepsize)
  proposed_v <- - leapfrog_result$v
  proposed_x <- leapfrog_result$x
  proposed_pdf <- pb$target(proposed_x)
  accept_ratio <- proposed_pdf - current_pdf
  # the acceptance ratio also features the "kinetic energy" term of the extended target
  # accept_ratio <- accept_ratio + sum(current_v^2) / 2 - sum(proposed_v^2) / 2
  accept_ratio <- accept_ratio + fast_dmvnorm_chol_inverse(matrix(proposed_v, nrow=1), rep(0, dimension_), pb$sigma_chol) -
    fast_dmvnorm_chol_inverse(matrix(current_v, nrow = 1), rep(0, dimension_), pb$sigma_chol)
  accept <- FALSE
  if (is.finite(accept_ratio)){
    accept <- (log(runif(1)) < accept_ratio)
  }
  if (accept){
    chain_state <- proposed_x
    current_pdf <- proposed_pdf
    accept <- TRUE
  }
  return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
}

hmc_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, nsteps, stepsize, pb){
  dimension_ <- length(chain_state1)
  current_v <- fast_rmvnorm(1, rep(0, dimension_), pb$precision_pi)
  leapfrog_result1 <- leapfrog_c(chain_state1, current_v, pb$precision_pi, pb$Sigma_pi, nsteps, stepsize)
  leapfrog_result2 <- leapfrog_c(chain_state2, current_v, pb$precision_pi, pb$Sigma_pi, nsteps, stepsize)
  proposed_v1 <- - leapfrog_result1$v
  proposed_x1 <- leapfrog_result1$x
  proposed_v2 <- - leapfrog_result2$v
  proposed_x2 <- leapfrog_result2$x
  # if we were a bit smarter, we would save logtarget(current_x) so as
  # to not re-evaluate it at every step
  proposed_pdf1 <- pb$target(proposed_x1)
  proposed_pdf2 <- pb$target(proposed_x2)
  accept_ratio1 <- proposed_pdf1 - current_pdf1
  accept_ratio2 <- proposed_pdf2 - current_pdf2
  # the acceptance ratio also features the "kinetic energy" term of the extended target
  accept_ratio1 <- accept_ratio1 + fast_dmvnorm_chol_inverse(matrix(proposed_v1, nrow=1), rep(0, dimension_), pb$sigma_chol) -
    fast_dmvnorm_chol_inverse(matrix(current_v, nrow = 1), rep(0, dimension_), pb$sigma_chol)
  accept_ratio2 <- accept_ratio2 + fast_dmvnorm_chol_inverse(matrix(proposed_v2, nrow=1), rep(0, dimension_), pb$sigma_chol) -
    fast_dmvnorm_chol_inverse(matrix(current_v, nrow = 1), rep(0, dimension_), pb$sigma_chol)
  accept <- 0
  logu <- log(runif(1)) # shared by both chains
  if (logu < accept_ratio1){
    accept <- accept + 1
    chain_state1 <- proposed_x1
    current_pdf1 <- proposed_pdf1
  } else {
    #
  }
  if (logu < accept_ratio2){
    accept <- accept + 1
    chain_state2 <- proposed_x2
    current_pdf2 <- proposed_pdf2
  } else {
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, current_pdf1 = current_pdf1, current_pdf2 = current_pdf2,
              accept = accept))
}

mh_single_kernel <- function(chain_state, current_pdf, stepsize = 1, pb){
  proposal_value <- chain_state + stepsize * rnorm(length(chain_state), mean = 0, sd = 1)
  proposal_pdf <- pb$target(proposal_value)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
  }
}

# max coupling of Gaussians, centered at mu1 and mu2, and with covariance matrix stepsize^2 * Identity
gaussian_max_coupling <- function(mu1, mu2, stddev){
  d <- length(mu1)
  x <- rnorm(d, mu1, stddev)
  if (sum(dnorm(x, mu1, stddev, log = TRUE)) + log(runif(1)) < sum(dnorm(x, mu2, stddev, log = TRUE))){
    return(cbind(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rnorm(d, mu2, stddev)
      reject <- (sum(dnorm(y, mu2, stddev, log = TRUE)) + log(runif(1)) < sum(dnorm(y, mu1, stddev, log = TRUE)))
    }
    return(cbind(x,y))
  }
}
# coupled MH kernel
mh_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, stepsize, pb){
  proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, stepsize)
  proposal1 <- proposal_value[,1]
  proposal2 <- proposal_value[,2]
  proposal_pdf1 <- pb$target(proposal1)
  proposal_pdf2 <- pb$target(proposal2)
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
    current_pdf1 <- proposal_pdf1
  }
  if (accept2){
    chain_state2 <- proposal2
    current_pdf2 <- proposal_pdf2
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, current_pdf1 = current_pdf1, current_pdf2 = current_pdf2))
}

# scaling prescribed in Roberts et al.
scaling_stepsize <- function(constant, dimension){
  return(constant * (1 / dimension)^(1/4))
}
#


samplemeetingtime <- function(omega, mh_stepsize, hmcconstant, dimension, max_iterations = Inf, pb){
  # HMC parameters
  hmc_stepsize <- scaling_stepsize(hmcconstant, dimension)
  # integration time is set to one
  hmc_nsteps <- 1 + floor(1 / hmc_stepsize)
  # Mixture kernels
  mixture_single_kernel <- function(chain_state, current_pdf){
    if (runif(1) < omega){
      return(mh_single_kernel(chain_state, current_pdf, mh_stepsize, pb))
    } else {
      return(hmc_single_kernel(chain_state, current_pdf, hmc_nsteps, hmc_stepsize, pb))
    }
  }
  #
  mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    if (runif(1) < omega){
      res_ <- mh_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, mh_stepsize, pb)
      res_$hmc_attempt <- 0
      return(res_)
    } else {
      res_ <- hmc_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, hmc_nsteps, hmc_stepsize, pb)
      res_$hmc_attempt <- 1
      return(res_)
    }
  }
  #
  chain_state1 <- pb$rinit(dimension)
  chain_state2 <- pb$rinit(dimension)
  current_pdf1 <- pb$target(chain_state1)
  current_pdf2 <- pb$target(chain_state2)
  sres1 <- mixture_single_kernel(chain_state1, current_pdf1)
  chain_state1 <- sres1$chain_state
  current_pdf1 <- sres1$current_pdf
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  hmc_accept <- 0
  hmc_attempt <- 0

  while (!finished && iter < max_iterations){
    iter <- iter + 1
    res_coupled_kernel <- mixture_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
    if (res_coupled_kernel$hmc_attempt == 1){
      hmc_attempt <- hmc_attempt + 1
      hmc_accept <- hmc_accept + res_coupled_kernel$accept
    }
    chain_state1 <- res_coupled_kernel$chain_state1
    chain_state2 <- res_coupled_kernel$chain_state2
    current_pdf1 <- res_coupled_kernel$current_pdf1
    current_pdf2 <- res_coupled_kernel$current_pdf2
    if (all(chain_state1 == chain_state2) && !meet){
      # recording meeting time tau
      meet <- TRUE
      meetingtime <- iter
    }
    # stop after tau steps
    if (iter >= meetingtime){
      finished <- TRUE
    }
  }
  return(list(meetingtime = meetingtime, iteration = iter, finished = finished, hmc_acceptrate = hmc_accept/(2*hmc_attempt)))
}

# dimensions <- c(10, 50, 100)
dimensions <- c(10, 50, 100, 200, 300)
nsamples <- 1000
df <- data.frame()
precisions <- list()
for (idim in seq_along(dimensions)){
  dimension <- dimensions[idim]
  omega <- 0.05
  hmcconstant <- .1
  mh_stepsize <- 1e-4
  hmc_stepsize <- scaling_stepsize(hmcconstant, dimension)
  hmc_nsteps <- 1 + floor(1 / hmc_stepsize)
  print(dimension)
  # sample meeting times
  results <- foreach(irep = 1:nsamples) %dorng% {
    precision_pi <- rWishart(1, dimension, diag(dimension))[,,1]
    # precisions[[idim]] <- precision_pi
    Sigma_pi <- solve(precision_pi)
    precision_chol <- t(chol(precision_pi))
    sigma_chol <- t(chol(Sigma_pi))
    mean_target <- rep(0, dimension)
    # log-density of multivariate Normal
    target <- function(x){
      return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), mean_target, precision_chol))
    }
    # gradient of log-density of multivariate Normal
    gradtarget <- function(x) gradlogtarget_c(x, mean_target, precision_pi)
    ## Starting from stationarity
    rinit <- function(dimension){
      return(fast_rmvnorm(1, mean_target, Sigma_pi))
    }
    pb <- list(rinit = rinit, target = target, gradtarget = gradtarget, Sigma_pi = Sigma_pi, precision_chol = precision_chol,
               sigma_chol = sigma_chol, mean_target = mean_target, precision_pi = precision_pi)
    samplemeetingtime(omega, mh_stepsize, hmcconstant, dimension, max_iterations = 1e4, pb)
  }
  meetingtimes <- sapply(results, function(x) x$meetingtime)
  print(summary(meetingtimes))
  df <- rbind(df, data.frame(idim = idim, dimension = dimension, hmc_nsteps = hmc_nsteps, hmc_stepsize = hmc_stepsize,
             irep = 1:nsamples, meetingtimes = meetingtimes, init_type = "target"))
  ## Starting outside of stationarity
  ## in that case, the initialization is from Normal(1,I)
  ## and then, 10 times in a row, the chain is propagated through the leap frog integrator
  ## without doing the MH acceptance step as in the usual HMC kernel

  # sample meeting times
  results <- foreach(irep = 1:nsamples) %dorng% {
    precision_pi <- rWishart(1, dimension, diag(dimension))[,,1]
    # precisions[[idim]] <- precision_pi
    Sigma_pi <- solve(precision_pi)
    precision_chol <- t(chol(precision_pi))
    sigma_chol <- t(chol(Sigma_pi))
    mean_target <- rep(0, dimension)
    # log-density of multivariate Normal
    target <- function(x){
      return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), mean_target, precision_chol))
    }
    # gradient of log-density of multivariate Normal
    gradtarget <- function(x) gradlogtarget_c(x, mean_target, precision_pi)
    rinit <- function(dimension){
      chain_state <- rnorm(dimension, mean = 1)
      for (iter in 1:10){
        current_v <- fast_rmvnorm(1, rep(0, dimension), precision_pi)
        leapfrog_result <- leapfrog_c(chain_state, current_v, precision_pi, Sigma_pi, hmc_nsteps, hmc_stepsize)
        chain_state <- leapfrog_result$x
      }
      return(chain_state)
    }
    pb <- list(rinit = rinit, target = target, gradtarget = gradtarget, Sigma_pi = Sigma_pi, precision_chol = precision_chol,
               sigma_chol = sigma_chol, mean_target = mean_target, precision_pi = precision_pi)
    samplemeetingtime(omega, mh_stepsize, hmcconstant, dimension, max_iterations = 1e4, pb)
  }
  meetingtimes <- sapply(results, function(x) x$meetingtime)
  df <- rbind(df, data.frame(idim = idim, dimension = dimension, hmc_nsteps = hmc_nsteps, hmc_stepsize = hmc_stepsize,
                             irep = 1:nsamples, meetingtimes = meetingtimes, init_type = "offset"))
  save(df, file = "scalingdimension.wishart.hmc.meetings.RData")
}

load("scalingdimension.wishart.hmc.meetings.RData")

# head(df)
# max(df$meetingtimes)
# library(dplyr)
# library(tidyr)
# df %>% filter(irep == 1)
#
# df.summary <- df %>% group_by(dimension, init_type) %>% summarise(mean_time = mean(meetingtimes))
#
# g <- ggplot(df.summary, aes(x = dimension, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.summary$dimension))) + xlab("dimension")
# g

