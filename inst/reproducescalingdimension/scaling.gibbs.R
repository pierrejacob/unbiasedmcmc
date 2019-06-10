# load packages
library(debiasedmcmc)
library(doParallel)
library(doRNG)
# library(ggthemes)
# setmytheme()
rm(list = ls())
set.seed(1)
registerDoParallel(cores = detectCores()-2)
#

get_problem <- function(dimension, iterate_d,
                        target_type='sparse',
                        init_type='target',
                        proposal_type='scaled'
){
  iterate_per_component <- iterate_d(dimension)

  if(target_type=='sparse'){
    Sigma_pi <- diag(1, dimension, dimension)
    alpha <- 0.5
    for (i in 1:dimension){
      for (j in 1:dimension){
        Sigma_pi[i,j] <- alpha^(abs(i-j))
      }
    }
  } else if (target_type=='dense'){
    Sigma_pi <- solve(rWishart(1, dimension, diag(dimension))[,,1])
  }

  Sigma_chol <- chol(Sigma_pi)
  Sigma_chol_inv <- solve(chol(Sigma_pi))
  m <- rep(0, dimension)
  target <- function(x) fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), m, Sigma_chol_inv)

  if (proposal_type == 'unscaled'){
    proposal_sd_vec <- rep(1, dimension)
  } else if (proposal_type == 'scaled') {
    proposal_sd_vec <- sqrt(diag(Sigma_pi))
  }

  ## Single Gibbs update and kernel
  gibbs_update <- function(component, proposal_sd){
    single_step_ <- function(chain_state, proposal_sd){
      increment <- rnorm(1, sd=proposal_sd)
      proposal_value <- chain_state
      proposal_value[component] <- proposal_value[component] + increment
      proposal_pdf <- target(proposal_value)
      current_pdf <- target(chain_state)
      accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
      if (accept){
        return(proposal_value)
      } else {
        return(chain_state)
      }
    }
    return(single_step_)
  }

  single_kernel <- function(chain_state){
    for (component in 1:dimension){
      proposal_sd <- proposal_sd_vec[[component]]
      step_ <- gibbs_update(component, proposal_sd)
      for (iterate in 1:iterate_per_component){
        chain_state <- step_(chain_state, proposal_sd)
      }
    }
    return(chain_state)
  }

  ## Coupled Gibbs update and kernel
  coupled_gibbs_update <- function(component, proposal_sd){
    coupled_kernel_ <- function(chain_state1, chain_state2){
      proposal_value <- rnorm_max_coupling(chain_state1[component],
                                           chain_state2[component],
                                           proposal_sd, proposal_sd)
      proposal1 <- chain_state1
      proposal1[component] <- proposal_value$xy[1]
      proposal2 <- chain_state2
      proposal2[component] <- proposal_value$xy[2]
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
    return(coupled_kernel_)
  }

  coupled_kernel <- function(chain_state1, chain_state2){
    for (component in 1:dimension){
      proposal_sd <- proposal_sd_vec[[component]]
      step_ <- coupled_gibbs_update(component, proposal_sd)
      for (iterate in 1:iterate_per_component){
        chain_states <- step_(chain_state1, chain_state2)
        chain_state1 <- chain_states$chain_state1
        chain_state2 <- chain_states$chain_state2
      }
    }
    identical_ <- all(chain_state1 == chain_state2)
    return(list(state1 = chain_state1, state2 = chain_state2, identical = identical_))
  }

  if (init_type == 'target'){
    rinit <- function() fast_rmvnorm(1, rep(0, dimension), Sigma_pi)
  } else if (init_type == 'offset'){
    rinit <- function() fast_rmvnorm(1, rep(1, dimension), diag(dimension))
  }

  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit, Sigma_pi = Sigma_pi))
}

## Sparse case
nsamples <- 1000
init_type_vec <- c('target','offset')
proposal_type_vec <- 'scaled'
run_parameters <- expand.grid(d=c(1, 50, 100, 200, 300),
                              iterates = 1,
                              target_type = 'sparse',
                              init_type=init_type_vec,
                              proposal_type=proposal_type_vec,
                              stringsAsFactors = FALSE)

# run_parameters

df <- data.frame()
for (i in 1:nrow(run_parameters)){
  d = run_parameters[i,'d']

  iterate = run_parameters[i,'iterates']
  iterate_d <- function(dimension) iterate

  target_type = run_parameters[i,'target_type']
  init_type = run_parameters[i,'init_type']
  proposal_type = run_parameters[i,'proposal_type']


  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    problem <- get_problem(d, iterate_d, target_type, init_type, proposal_type)
    sample_meetingtime(problem$single_kernel, problem$coupled_kernel, problem$rinit)
  }
  meetingtime <- sapply(c_chains_, function(x) x$meetingtime)

  cat("dim:", d, ", iterate:", iterate, ", target:", target_type, ", init:", init_type, "\n")
  print(summary(meetingtime))
  q5 <- as.numeric(quantile(meetingtime, probs = 0.05))
  q95 <- as.numeric(quantile(meetingtime, probs = 0.95))
  df <- rbind(df, data.frame(d=d, iterate=iterate,
                             target_type=target_type,
                             init_type=init_type,
                             proposal_type=proposal_type,
                             mean_time=mean(meetingtime),
                             q5 = q5, q95 = q95,
                             max_time=max(meetingtime),
                             median_time=median(meetingtime)))
  save(nsamples, run_parameters, df, file = "scalingdimension.gibbs.sparse.RData")
}
load(file = "scalingdimension.gibbs.sparse.RData")
df.sparse <- df

# Now, the dense case
nsamples <- 1000
run_parameters <- expand.grid(d=floor(seq(from = 1, to = 9, length.out = 5)),
                              iterates = 1,
                              target_type = 'dense',
                              proposal_type=proposal_type_vec,
                              stringsAsFactors = FALSE)

run_parameters
df <- data.frame()
for (i in 1:nrow(run_parameters)){
  d = run_parameters[i,'d']

  iterate = run_parameters[i,'iterates']
  iterate_d <- function(dimension) iterate

  target_type = run_parameters[i,'target_type']
  proposal_type = run_parameters[i,'proposal_type']

  init_type = 'target'
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    problem <- get_problem(d, iterate_d, target_type, 'target', proposal_type)
    problem$rinit <- function() fast_rmvnorm(1, rep(0, d), problem$Sigma_pi)
    sample_meetingtime(problem$single_kernel, problem$coupled_kernel, problem$rinit, max_iterations = 5e5)
  }
  meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
  ninf <- sum(is.infinite(meetingtime))
  meetingtime <- meetingtime[is.finite(meetingtime)]
  q5 <- as.numeric(quantile(meetingtime, probs = 0.05))
  q95 <- as.numeric(quantile(meetingtime, probs = 0.95))
  df.init_target <- data.frame(d=d, iterate=iterate,
                               target_type=target_type,
                               init_type=init_type,
                               proposal_type=proposal_type,
                               mean_time=mean(meetingtime),
                               q5 = q5, q95 = q95,
                               ninf = ninf,
                               max_time=max(meetingtime),
                               median_time=median(meetingtime))
  cat("dim:", d, ", iterate:", iterate, ", target:", target_type, ", init:", init_type, "\n")
  print(summary(meetingtime))
  init_type = 'offset'
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    problem <- get_problem(d, iterate_d, target_type, 'target', proposal_type)
    problem$rinit <- function() fast_rmvnorm(1, rep(1, d), diag(d))
    sample_meetingtime(problem$single_kernel, problem$coupled_kernel, problem$rinit, max_iterations = 5e5)
  }
  meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
  ninf <- sum(is.infinite(meetingtime))
  meetingtime <- meetingtime[is.finite(meetingtime)]
  q5 <- as.numeric(quantile(meetingtime, probs = 0.05))
  q95 <- as.numeric(quantile(meetingtime, probs = 0.95))
  df.init_offset <- data.frame(d=d, iterate=iterate,
                               target_type=target_type,
                               init_type=init_type,
                               proposal_type=proposal_type,
                               mean_time=mean(meetingtime),
                               ninf = ninf,
                               q5 = q5, q95 = q95,
                               max_time=max(meetingtime),
                               median_time=median(meetingtime))

  cat("dim:", d, ", iterate:", iterate, ", target:", target_type, ", init:", init_type, "\n")
  print(summary(meetingtime))
  df <- rbind(df, df.init_target, df.init_offset)
  save(nsamples, run_parameters, df, file = "scalingdimension.gibbs.dense.RData")
}
load("scalingdimension.gibbs.dense.RData")
df.dense <- df
df.dense
#
# library(dplyr)
# g <- ggplot(df.sparse, aes(x = d, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.sparse$d))) + xlab("dimension")
# g
#
#
# # ggsave(filename = "scalingdimension.gibbs.sparse.pdf", plot = g, width = 8, height = 6)
#
# g <- ggplot(df.dense, aes(x = d, y = median_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.dense$d))) + xlab("dimension")
# g <- g + scale_y_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5))
# g
# # ggsave(filename = "scalingdimension.gibbs.dense.pdf", plot = g, width = 8, height = 6)
#
# g <- ggplot(df.dense, aes(x = d, y = median_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("median meeting time")
# g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.dense$d))) + xlab("dimension")
# g

# g <- ggplot(df.dense, aes(x = d, ymin = q5, ymax = q95, group = init_type, fill = init_type)) + geom_ribbon(alpha = 0.5)
#
#             , linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique(df.dense$d))) + xlab("dimension")
# g

