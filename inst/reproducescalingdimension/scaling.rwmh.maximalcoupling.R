# load packages
library(debiasedmcmc)
library(doParallel)
library(doRNG)
library(dplyr)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores()-2)
#
library(ggthemes)


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

  # Markov kernel of the chain
  single_kernel <- function(chain_state){
    proposal_value <- fast_rmvnorm(1, chain_state, Sigma_proposal)
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
  coupled_kernel <- function(chain_state1, chain_state2){
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

  if (init_type=='target'){
    rinit <- function() fast_rmvnorm(1, rep(0, dimension), Sigma_pi)
  } else if (init_type=='offset'){
    rinit <- function() fast_rmvnorm(1, rep(1, dimension), diag(dimension))
  }
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit))
}

nsamples <- 1000

dim_vec <- floor(seq(from = 1, to = 9, length.out = 5))
target_type_vec <- c('sparse','dense')[2]
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

  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    problem <- get_problem(d, target_type, init_type, scale_type=scale_type)
    coupled_chains(problem$single_kernel, problem$coupled_kernel, problem$rinit)
  }
  meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
  cat("dim ", d, ", target ", target_type, ", init ", init_type, ", scale ", scale_type, "\n")
  print(summary(meetingtime))
  q5 <- as.numeric(quantile(meetingtime, probs = 0.05))
  q95 <- as.numeric(quantile(meetingtime, probs = 0.95))
  df <- rbind(df, data.frame(d = d,
                             target_type = target_type,
                             init_type = init_type,
                             scale_type = scale_type,
                             mean_time = mean(meetingtime), max_time = max(meetingtime),
                             q5 = q5, q95 = q95,
                             median_time = median(meetingtime)))
  save(nsamples, run_parameters, df, file = "scalingdimension.rwmh.maxcoupling.RData")
}

load(file = "scalingdimension.rwmh.maxcoupling.RData")
df

# g <- ggplot(df %>% filter(target_type == "sparse"), aes(x = d, y = mean_time, group = init_type, colour = factor(init_type))) + geom_line() + ylab("average meeting time")
# g <- g + scale_color_colorblind() + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "sparse"))$d)))
# g + scale_y_log10()

g <- ggplot(df %>% filter(target_type == "dense"), aes(x = d, y = mean_time, group = init_type, linetype = factor(init_type))) + geom_line() + ylab("average meeting time")
g <- g + scale_linetype("initialization:") + scale_x_continuous(breaks = sort(unique((df %>% filter(target_type == "dense"))$d)))
g <- g + scale_y_log10(breaks = c(10,100,1e3, 1e4)) + xlab("dimension")
g
# ggsave(filename = "scalingdimension.rwmh.maxcoupling.pdf", plot = g, width = 8, height = 6)
#
