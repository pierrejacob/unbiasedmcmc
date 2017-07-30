# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
setwd("~/Dropbox/PolyaGammaResults/dimension/")

get_problem <- function(dimension, iterated_d){
  iterate_per_component <- iterated_d(dimension)
  Sigma_pi <- diag(1, dimension, dimension)
  alpha <- 0.5
  for (i in 1:dimension){
    for (j in 1:dimension){
      Sigma_pi[i,j] <- alpha^(abs(i-j))
    }
  }

  target <- function(x) fast_dmvnorm(matrix(x, nrow = 1), rep(0, dimension), Sigma_pi)
  ## now try the MH within Gibbs way
  gibbs_update <- function(component){
    single_step_ <- function(chain_state){
      increment <- rnorm(1)
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
      step_ <- gibbs_update(component)
      for (iterate in 1:iterate_per_component){
        chain_state <- step_(chain_state)
      }
    }
    return(chain_state)
  }
  coupled_gibbs_update <- function(component){
    coupled_kernel_ <- function(chain_state1, chain_state2){
      proposal_value <- gaussian_max_coupling(chain_state1[component], chain_state2[component], diag(1,1,1), diag(1,1,1))
      proposal1 <- chain_state1
      proposal1[component] <- proposal_value[,1]
      proposal2 <- chain_state2
      proposal2[component] <- proposal_value[,2]
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
      step_ <- coupled_gibbs_update(component)
      for (iterate in 1:iterate_per_component){
        chain_states <- step_(chain_state1, chain_state2)
        chain_state1 <- chain_states$chain_state1
        chain_state2 <- chain_states$chain_state2
      }
    }
    return(chain_states)
  }
  rinit <- function() fast_rmvnorm(1, rep(0, dimension), Sigma_pi)
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel, rinit = rinit))
}

iterate_d <- function(dimension) floor(1 + log(dimension))

# df <- data.frame()
# for (d in floor(seq(from = 1, to = 100, length.out = 30))){
#   problem <- get_problem(d, iterate_d)
#   nsamples <- 1000
#   c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#     coupled_chains(problem$single_kernel, problem$coupled_kernel, problem$rinit, max_iterations = 1e4)
#   }
#   meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
#   cat("dimension ", d, "\n")
#   print(summary(meetingtime))
#   df <- rbind(df, data.frame(d = d, mean_time = mean(meetingtime), max_time = max(meetingtime),
#                              median_time = median(meetingtime)))
#   save(df, file = "dimension.scaling.gibbs.RData")
# }
load(file = "dimension.scaling.gibbs.RData")
g <- ggplot(df, aes(x = d, y = mean_time)) + geom_line() + scale_x_continuous(breaks = c(1,25, 50,75, 100)) + ylab("average meeting time")
g <- g + scale_y_continuous(breaks = c(1, 10, 20), limits = c(0,20))
g
ggsave(filename = "dimension.scaling.gibbs.pdf", plot = g, width = 7, height = 7)


iterate_d <- function(dimension) floor(1 + dimension)
df <- data.frame()
for (d in floor(seq(from = 1, to = 100, length.out = 30))){
  problem <- get_problem(d, iterate_d)
  nsamples <- 1000
  c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
    coupled_chains(problem$single_kernel, problem$coupled_kernel, problem$rinit, max_iterations = 1e4)
  }
  meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
  cat("dimension ", d, "\n")
  print(summary(meetingtime))
  df <- rbind(df, data.frame(d = d, mean_time = mean(meetingtime), max_time = max(meetingtime),
                             median_time = median(meetingtime)))
  save(df, file = "dimension.scaling.gibbs.iterlinear.RData")
}
load(file = "dimension.scaling.gibbs.iterlinear.RData")
g <- ggplot(df, aes(x = d, y = mean_time)) + geom_line() + scale_x_continuous(breaks = c(1,25, 50,75, 100)) + ylab("average meeting time")
g <- g + scale_y_continuous(breaks = c(1, 10, 20), limits = c(0,50))
g
# ggsave(filename = "dimension.scaling.gibbs.pdf", plot = g, width = 7, height = 7)

