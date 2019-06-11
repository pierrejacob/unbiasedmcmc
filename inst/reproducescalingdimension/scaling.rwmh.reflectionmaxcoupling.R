# load packages
library(unbiasedmcmc)
library(doParallel)
library(doRNG)
# setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores()-2)
#

## obtain functions to generate MH chains targeting a multivariate Normal density
get_problem <- function(dimension,
                        target_type='sparse',
                        init_type='target',
                        scale_type='unscaled'){
  mean_pi <- rep(0, dimension)
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
  Sigma_pi_chol <- chol(Sigma_pi)
  Sigma_pi_chol_inv <- solve(Sigma_pi_chol)

  target <- function(x) fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), mean_pi, Sigma_pi_chol_inv)

  # Set proposal distribution scaling
  if(scale_type=='unscaled'){
    Sigma_proposal <- Sigma_pi
  } else if(scale_type=='optimal'){
    Sigma_proposal <- Sigma_pi/dimension
  }
  Sigma_proposal_chol <- chol(Sigma_proposal)
  Sigma_proposal_chol_inv <- solve(Sigma_proposal_chol)


  # Markov kernel, taking a state and associated log-pdf of the target distribution
  single_kernel <- function(state){
    chain_state <- state$chain_state
    current_pdf <- state$current_pdf
    proposal_value <- fast_rmvnorm_chol(1, chain_state, Sigma_proposal_chol)
    proposal_pdf <- target(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
  }

  # Markov kernel of the coupled chain
  coupled_kernel <- function(state1, state2){
    chain_state1 <- state1$chain_state; current_pdf1 <- state1$current_pdf
    chain_state2 <- state2$chain_state; current_pdf2 <- state2$current_pdf
    proposal_value <- rmvnorm_reflectionmax(chain_state1, chain_state2, Sigma_proposal_chol, Sigma_proposal_chol_inv)
    proposal1 <- proposal_value$xy[,1]; proposal_pdf1 <- target(proposal1)
    if (proposal_value$identical){
      proposal2 <- proposal1; proposal_pdf2 <- proposal_pdf1
    } else {
      proposal2 <- proposal_value$xy[,2]; proposal_pdf2 <- target(proposal2)
    }
    logu <- log(runif(1))
    accept1 <- FALSE; accept2 <- FALSE
    if (is.finite(proposal_pdf1)) accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    if (is.finite(proposal_pdf2)) accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    if (accept1){
      chain_state1 <- proposal1
      current_pdf1 <- proposal_pdf1
      }
    if (accept2){
      chain_state2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    identical_ <- proposal_value$identical && accept1 && accept2
    return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
                state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
                identical = identical_))
  }
  if (init_type=='target'){
    rinit <- function(){
      state_ <- fast_rmvnorm_chol(1, mean_pi, Sigma_pi_chol)[1,]
      pdf_ <- target(state_)
      return(list(chain_state = state_, current_pdf = pdf_))
    }
  } else if (init_type=='offset'){
    rinit <- function(){
      state_ <- fast_rmvnorm_chol(1, rep(1, dimension), diag(dimension))[1,]
      pdf_ <- target(state_)
      return(list(chain_state = state_, current_pdf = pdf_))
    }
  }
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit,
              target = target, Sigma_proposal = Sigma_proposal))
}
#

nsamples <- 1000
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
    sample_meetingtime(problem$single_kernel, problem$coupled_kernel, problem$rinit)
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

# library(dplyr)
# g <- ggplot(df %>% filter(target_type == "dense"), aes(x = d, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g +  scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
# g <- g + xlab("dimension")
# g
# ggsave(filename = "scalingdimension.rwmh.reflmaxcoupling.pdf", plot = g, width = 8, height = 6)
# g
