library(debiasedmcmc)
rm(list = ls())
set.seed(18)
setmytheme()
registerDoParallel(cores = detectCores())

# number of leapfrog steps
nsteps <- 20
#

cppFunction('
            Eigen::VectorXd gradlogtarget_c(const Eigen::VectorXd & x, const Eigen::VectorXd & mean, const Eigen::MatrixXd & precision){
return(- precision * (x - mean));
            }
            ', depends = "RcppEigen")

get_mvnormal <- function(dimension, mean_target, Sigma_target){
  # compute precision matrix
  precision <- solve(Sigma_target)
  # compute Cholesky factor of precision matrix
  precision_chol <- t(chol(precision))
  # log-density of multivariate Normal
  logtarget <- function(x){
    return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), mean_target, precision_chol))
  }
  # gradient of log-density of multivariate Normal
  # gradlogtarget <- function(x){
  #   return((- precision %*% matrix(x - mean_target, ncol = 1))[,1])
  # }
  gradlogtarget <- function(x) gradlogtarget_c(x, mean_target, precision)
  return(list(logtarget = logtarget, gradlogtarget = gradlogtarget))
}

get_hmc_kernel <- function(logtarget, gradlogtarget, stepsize, nsteps, dimension){
  # leap frog integrator
  # note that some papers use the notation U for - logtarget, so that there are minus signs everywhere
  leapfrog <- function(x, v){
    # xtraj <- matrix(nrow = nsteps+1, ncol = length(x))
    # xtraj[1,] <- x
    v <- v + stepsize * gradlogtarget(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * v
      # xtraj[step+1,] <- x
      if (step != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    # we could negate the momentum but we don't use it here
    # return(list(x = x, v = v, xtraj = xtraj))
    return(list(x = x, v = v))
  }
  # One step of HMC
  kernel <- function(chain_state, iteration){
    current_v <- rnorm(dimension) # velocity or momentum
    leapfrog_result <- leapfrog(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x
    # if we were a bit smarter, we would save logtarget(chain_state) so as
    # to not re-evaluate it at every step
    accept_ratio <- logtarget(proposed_x) - logtarget(chain_state)
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio <- accept_ratio + sum(current_v^2) / 2 - sum(proposed_v^2) / 2
    # we store the entire x-trajectory here, for plotting purposes
    # the first row stores the current value (which is also the last value of the latest accepted trajectory)
    xtraj <- matrix(nrow = nsteps+1, ncol = dimension)
    accept <- FALSE
    if (log(runif(1)) < accept_ratio){
      chain_state <- proposed_x
      current_v <- proposed_v
      xtraj <- leapfrog_result$xtraj
      accept <- TRUE
    } else {
      # if we reject, the entire x-trajectory is boring
      for (istep in 1:(nsteps+1)){
        xtraj[istep,] <- chain_state
      }
    }
    return(list(chain_state = chain_state, xtraj = xtraj, accept = accept))
  }
  # One step of HMC
  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    current_v <- rnorm(dimension) # velocity or momentum, shared by both chains
    leapfrog_result1 <- leapfrog(chain_state1, current_v)
    leapfrog_result2 <- leapfrog(chain_state2, current_v)
    proposed_v1 <- - leapfrog_result1$v
    proposed_x1 <- leapfrog_result1$x
    proposed_v2 <- - leapfrog_result2$v
    proposed_x2 <- leapfrog_result2$x
    # if we were a bit smarter, we would save logtarget(current_x) so as
    # to not re-evaluate it at every step
    accept_ratio1 <- logtarget(proposed_x1) - logtarget(chain_state1)
    accept_ratio2 <- logtarget(proposed_x2) - logtarget(chain_state2)
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio1 <- accept_ratio1 + sum(current_v^2) / 2 - sum(proposed_v1^2) / 2
    accept_ratio2 <- accept_ratio2 + sum(current_v^2) / 2 - sum(proposed_v2^2) / 2
    # we store the entire x-trajectory here, for plotting purposes
    # the first row stores the current value (which is also the last value of the latest accepted trajectory)
    xtraj1 <- matrix(nrow = nsteps+1, ncol = dimension)
    xtraj2 <- matrix(nrow = nsteps+1, ncol = dimension)
    logu <- log(runif(1)) # shared by both chains
    if (logu < accept_ratio1){
      chain_state1 <- proposed_x1
      current_v1 <- proposed_v1
      xtraj1 <- leapfrog_result1$xtraj
    } else {
      # if we reject, the entire x-trajectory is boring
      for (istep in 1:(nsteps+1)){
        xtraj1[istep,] <- chain_state1
      }
    }
    if (logu < accept_ratio2){
      chain_state2 <- proposed_x2
      current_v2 <- proposed_v2
      xtraj2 <- leapfrog_result2$xtraj
    } else {
      for (istep in 1:(nsteps+1)){
        xtraj2[istep,] <- chain_state2
      }
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, xtraj1 = xtraj1, xtraj2 = xtraj2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}

get_mh_kernel <- function(logtarget, Sigma_proposal, dimension){
  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)
  # single kernel
  kernel <- function(chain_state, iteration){
    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
    proposal_pdf <- logtarget(proposal_value)
    current_pdf <- logtarget(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value))
    } else {
      return(list(chain_state = chain_state))
    }
  }
  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    # print("max coupling")
    # print(Sigma2_chol[1:5,1:5])
    # print(Sigma2_chol_inv[1:5,1:5])
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    # print(chain_state1)
    # print(chain_state2)
    # print(Sigma_proposal[1:5,1:5])
    # print(sum((chain_state1 - chain_state2)^2))
    # proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    # print("done")
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)
    current_pdf1 <- logtarget(chain_state1)
    current_pdf2 <- logtarget(chain_state2)
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
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}

get_problem <- function(dimension, var_prop_mh, effectiveTime){
  stepsize <- effectiveTime / nsteps
  Sigma_pi <- diag(1, dimension, dimension)
  alpha <- 0.5
  for (i in 1:dimension){
    for (j in 1:dimension){
      Sigma_pi[i,j] <- alpha^(abs(i-j))
    }
  }
  Sigma_chol <- chol(Sigma_pi)
  # Sigma_chol_inv <- solve(chol(Sigma_pi))
  mean_target <- rep(0, dimension)
  target <- get_mvnormal(dimension, mean_target, Sigma_pi)
  hmc_kernel <- get_hmc_kernel(target$logtarget, target$gradlogtarget, stepsize, nsteps, dimension)
  rinit <- function() fast_rmvnorm_chol(1, mean_target, Sigma_chol)
  omega <- 0.1
  Sigma_proposal <- var_prop_mh * diag(1, dimension, dimension)
  mh_kernel <- get_mh_kernel(target$logtarget, Sigma_proposal, dimension)
  # Mixture kernels
  mixture_kernel <- function(chain_state, iter){
    if (runif(1) < omega){
      return(mh_kernel$kernel(chain_state, iter))
    } else {
      return(hmc_kernel$kernel(chain_state, iter))
    }
  }
  #
  mixture_coupled_kernel <- function(chain_state1, chain_state2, iter){
    if (runif(1) < omega){
      return(mh_kernel$coupled_kernel(chain_state1, chain_state2, iter))
    } else {
      return(hmc_kernel$coupled_kernel(chain_state1, chain_state2, iter))
    }
  }
  return(list(single_kernel = mixture_kernel, coupled_kernel = mixture_coupled_kernel, rinit = rinit,
              target = target, mean_target = mean_target, Sigma_target = Sigma_pi))
}

samplemeetingtime <- function(single_kernel, coupled_kernel, rinit, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  sres1 <- single_kernel(chain_state1, 0)
  chain_state1 <- sres1$chain_state
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      sres1 <- single_kernel(chain_state1, iter)
      chain_state1 <- sres1$chain_state
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, iter)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    # stop after max(m, tau) steps
    if (iter >= meetingtime){
      finished <- TRUE
    }
  }
  return(list(meetingtime = meetingtime, iteration = iter, finished = finished))
}

# effectiveTime <- 3*pi/4
effectiveTime <- pi/2
pb <- get_problem(50, 1e-5, effectiveTime)

# x <- pb$rinit()
# precision <- solve(pb$Sigma_target)
# summary(pb$target$gradlogtarget(x) - gradlogtarget_c(x, pb$mean_target, precision))

nrep <- 10
max_iterations <- 1e3
cc <-  samplemeetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit, max_iterations = max_iterations)
cc

# ues <- foreach(irep = 1:nrep) %dorng% {
#  samplemeetingtime(pb$single_kernel, pb$coupled_kernel, pb$rinit, max_iterations = max_iterations)
# }
#
# sapply(ues, function(x) x$meetingtime)

##
# dimensions <- c(1, 10, 50, 100)
dimensions <- c(1, 50, 100, 200, 300)
var_prop_mh <- c(1e-5)
effectiveTimes <- c(pi/4, pi/3, pi/2)
# effectiveTimes <- c(pi/2)
nsamples <- 100
df <- data.frame()
for (effectiveTime in effectiveTimes){
  for (d in dimensions){
    problem <- get_problem(d, var_prop_mh, effectiveTime)
    cc_ <-  foreach(irep = 1:nsamples) %dorng% {
      samplemeetingtime(problem$single_kernel, problem$coupled_kernel, problem$rinit, max_iterations = 1e4)
    }
    meetingtime <- sapply(cc_, function(x) x$meetingtime)
    cat("dimension ", d, "\n")
    print(summary(meetingtime))
    df <- rbind(df, data.frame(d = d, effectiveTime = effectiveTime, mean_time = mean(meetingtime), max_time = max(meetingtime),
                               median_time = median(meetingtime)))
    save(effectiveTimes, dimensions, df, file = "scalingdimension.hmc.RData")
  }
}
#
load(file = "scalingdimension.hmc.RData")
df


# #
library(ggthemes)
g <- ggplot(df, aes(x = d, y = mean_time, group = effectiveTime, colour = factor(effectiveTime))) + geom_line() + geom_point() + ylab("average meeting time")
g <- g + scale_color_colorblind() + scale_x_continuous(breaks = dimensions) # + theme(legend.position = "none")
g <- g + xlab("dimension")
g

