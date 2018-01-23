# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#


## try coupling a MH for multivariate Normal directly
get_problem <- function(dimension){
  Sigma_pi <- diag(1, dimension, dimension)
  alpha <- 0.5
  for (i in 1:dimension){
    for (j in 1:dimension){
      Sigma_pi[i,j] <- alpha^(abs(i-j))
    }
  }
  target <- function(x) fast_dmvnorm(matrix(x, nrow = 1), rep(0, dimension), Sigma_pi)
  #
  Sigma_proposal <- Sigma_pi
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
  rinit <- function() fast_rmvnorm(1, rep(0, dimension), Sigma_pi)
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit))
}

nsamples <- 100

df <- data.frame()
for (d in floor(seq(from = 1, to = 15, length.out = 15))){
  problem <- get_problem(d)
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(problem$single_kernel, problem$coupled_kernel, problem$rinit)
  }
  meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
  cat("dimension ", d, "\n")
  print(summary(meetingtime))
  df <- rbind(df, data.frame(d = d, mean_time = mean(meetingtime), max_time = max(meetingtime),
                             median_time = median(meetingtime)))
  save(df, file = "scalingdimension.rwmh.scaling2.RData")
}

# load(file = "dimension.scaling.rwmh.RData")
# g <- ggplot(df, aes(x = d, y = mean_time)) + geom_line() + scale_x_continuous(breaks = seq(from = 1, to = 15, by = 2)) + ylab("average meeting time")
# g <- g + scale_y_log10() + geom_point()
# g
# ggsave(filename = "dimension.scaling.rwmh1.pdf", plot = g, width = 5, height = 5)

## try with optimal scaling
get_problem <- function(dimension){
  Sigma_pi <- diag(1, dimension, dimension)
  alpha <- 0.5
  for (i in 1:dimension){
    for (j in 1:dimension){
      Sigma_pi[i,j] <- alpha^(abs(i-j))
    }
  }
  target <- function(x) fast_dmvnorm(matrix(x, nrow = 1), rep(0, dimension), Sigma_pi)
  # Sigma_proposal <- diag(1/dimension, dimension, dimension)
  Sigma_proposal <- Sigma_pi/dimension
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
  rinit <- function() fast_rmvnorm(1, rep(0, dimension), Sigma_pi)
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit))
}

df <- data.frame()
for (d in floor(seq(from = 1, to = 10, length.out = 10))){
  problem <- get_problem(d)
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(problem$single_kernel, problem$coupled_kernel, problem$rinit)
  }
  meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
  cat("dimension ", d, "\n")
  print(summary(meetingtime))
  df <- rbind(df, data.frame(d = d, mean_time = mean(meetingtime), max_time = max(meetingtime),
                             median_time = median(meetingtime)))
  save(df, file = "scalingdimension.rwmh.scaling1.RData")
}

# load(file = "dimension.scaling.rwmh.optimal.RData")
# df <- df %>% filter(is.finite(max_time))
# g <- ggplot(df, aes(x = d, y = mean_time)) + geom_line() + scale_x_continuous(breaks = unique(df$d)) + ylab("average meeting time")
# g <- g + scale_y_log10(breaks = c(10, 100, 1000, 10000)) + geom_point()
# g
# ggsave(filename = "dimension.scaling.rwmh2.pdf", plot = g, width = 5, height = 5)

