# load packages
library(debiasedmcmc)
library(doParallel)
library(doRNG)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores()-2)
#

## try coupling a MH for multivariate Normal directly
get_problem <- function(dimension,
                        target_type='sparse',
                        init_type='target',
                        scale_type='unscaled'){
  Sigma_pi <- diag(1, dimension, dimension)

  if(target_type=='sparse'){
    alpha <- 0.5
    for (i in 1:dimension){
      for (j in 1:dimension){
        Sigma_pi[i,j] <- alpha^(abs(i-j))
      }
    }
  } else if (target_type=='dense'){
    Sigma_pi <- solve(rWishart(1, dimension, diag(dimension))[,,1])
    # Note that in Hoffman and Gelman (2014), the NUTS paper from which
    # we take our Wishart example, they use the same A matrix in all
  }
  target <- function(x) fast_dmvnorm(matrix(x, nrow = 1), rep(0, dimension), Sigma_pi)

  # Set proposal distribution scaling
  if(scale_type=='unscaled'){
    Sigma_proposal <- Sigma_pi
  } else if(scale_type=='optimal'){
    Sigma_proposal <- Sigma_pi/dimension
  }
  Sigma_chol <- chol(Sigma_proposal)
  inv_Sigma_chol <- solve(Sigma_chol)


  # Markov kernel of the chain
  single_kernel <- function(chain_state, current_pdf){
    proposal_value <- fast_rmvnorm(1, chain_state, Sigma_proposal)
    proposal_pdf <- target(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
  }

  # Markov kernel of the coupled chain
  coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    proposal_value <- debiasedmcmc:::rnorm_reflectionmax_(chain_state1, chain_state2, Sigma_chol, inv_Sigma_chol)
    proposal1 <- proposal_value$xy[,1]
    proposal2 <- proposal_value$xy[,2]
    if (proposal_value$identical){
      proposal2 <- proposal1
    }
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- target(proposal2)
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

  if (init_type=='target'){
    rinit <- function() fast_rmvnorm(1, rep(0, dimension), Sigma_pi)[1,]
  } else if (init_type=='offset'){
    rinit <- function() fast_rmvnorm(1, rep(1, dimension), diag(dimension))[1,]
  }
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit,
              target = target, Sigma_proposal = Sigma_proposal))
}

unbiasedestimator_ <- function(target, single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 1, max_iterations = Inf){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  current_pdf1 <- target(chain_state1)
  current_pdf2 <- target(chain_state2)
  # mcmcestimator computes the sum of h(X_t) for t=k,...,m
  mcmcestimator <- h(chain_state1)
  dimh <- length(mcmcestimator)
  if (k > 0){
    mcmcestimator <- rep(0, dimh)
  }
  # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
  correction <- rep(0, dimh)
  sres1 <- single_kernel(chain_state1, current_pdf1)
  chain_state1 <- sres1$chain_state
  current_pdf1 <- sres1$current_pdf
  if (k == 0){
    correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
  }
  if (k <= 1 && m >= 1){
    mcmcestimator <- mcmcestimator + h(chain_state1)
  }
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      sres1 <- single_kernel(chain_state1, current_pdf1)
      chain_state1 <- sres1$chain_state
      current_pdf1 <- sres1$current_pdf
      chain_state2 <- chain_state1
      if (k <= iter && iter <= m){
        mcmcestimator <- mcmcestimator + h(chain_state1)
      }
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
      chain_state1 <- res_coupled_kernel$chain_state1
      current_pdf1 <- res_coupled_kernel$current_pdf1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf2 <- res_coupled_kernel$current_pdf2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
      if (k <= iter){
        if (iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
        }
        correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
      }
    }
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  mcmcestimator <- mcmcestimator / (m - k + 1)
  uestimator <- mcmcestimator + correction
  return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}


nsamples <- 100
dim_vec <- floor(seq(from = 1, to = 50, length.out = 5))
target_type_vec <- c('dense')
init_type_vec <- c('target','offset')
scale_type_vec <- c('optimal')
run_parameters <- expand.grid(d=dim_vec,
                              target_type=target_type_vec,
                              init_type=init_type_vec,
                              scale_type=scale_type_vec,
                              stringsAsFactors = FALSE)

df <- data.frame()
for (i in 1:nrow(run_parameters)){
  d = run_parameters[i,'d']
  target_type = run_parameters[i,'target_type']
  init_type = run_parameters[i,'init_type']
  scale_type = run_parameters[i,'scale_type']
  res_ <- foreach(irep = 1:nsamples) %dorng% {
    problem <- get_problem(d, target_type, init_type, scale_type=scale_type)
    unbiasedestimator_(problem$target, problem$single_kernel, problem$coupled_kernel, problem$rinit, max_iterations = Inf)
  }
  meetingtime <- sapply(res_, function(x) x$meetingtime)
  cat("dim ", d, ", target ", target_type, ", init ", init_type, ", scale ", scale_type, "\n")
  print(summary(meetingtime))
  df <- rbind(df, data.frame(d = d,
                             target_type = target_type,
                             init_type = init_type,
                             scale_type = scale_type,
                             mean_time = mean(meetingtime), max_time = max(meetingtime),
                             median_time = median(meetingtime)))
  save(nsamples, run_parameters, df, file = "scalingdimension.rwmh.reflmaxcoupling.RData")
}

load(file = "scalingdimension.rwmh.reflmaxcoupling.RData")
df

library(dplyr)
# g <- ggplot(df %>% filter(target_type == "sparse"), aes(x = d, y = mean_time, group = init_type, colour = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g + scale_color_colorblind() + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "sparse"))$d)))
# g + ylim(0,2000)

g <- ggplot(df %>% filter(target_type == "dense"), aes(x = d, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
g <- g +  scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
g <- g + xlab("dimension")
g
# ggsave(filename = "scalingdimension.rwmh.reflmaxcoupling.pdf", plot = g, width = 8, height = 6)
# g

#
# dimension <- 5
# pb <- get_problem(dimension, target_type = "dense", init_type ='offset')
# unbiasedestimator_(pb$target, pb$single_kernel, pb$coupled_kernel, pb$rinit, max_iterations = 1e4)$meetingtime
#
# nrep <- 100
# dim_vec <- floor(seq(from = 5, to = 50, length.out = 5))
# df <- data.frame()
# for (dimension in dim_vec){
#   print(dimension)
#   pb <- get_problem(dimension, target_type = "sparse", scale_type = "optimal", init_type = "offset")
#   res_ <- foreach(irep = 1:nrep) %dorng% {
#     unbiasedestimator_(pb$target, pb$single_kernel, pb$coupled_kernel, pb$rinit, max_iterations = Inf)
#   }
#   meetings <- sapply(res_, function(x) x$meetingtime)
#   df <- rbind(df, data.frame(d = dimension,
#                              mean_time = mean(meetings),
#                              max_time = max(meetings),
#                              q90 = as.numeric(quantile(meetings, probs = 0.9)),
#                              median_time = median(meetings)))
#
# }
# #
# head(df)
# #
# ggplot(df, aes(x = d, y = mean_time)) + geom_line() + geom_point()
# ggplot(df, aes(x = d, y = q90, group = coupling, colour = coupling)) + geom_line() + viridis::scale_color_viridis(discrete=T)

