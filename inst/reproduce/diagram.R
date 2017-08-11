# load packages
library(debiasedmcmc)
library(ggthemes)
library(lars)
library(latex2exp)
setmytheme()
rm(list = ls())
set.seed(81)


## target distribution
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

sd_proposal <- 3
# Markov kernel of the chain
single_kernel <- function(chain_state){
  proposal_value <- rnorm(1, mean=chain_state, sd=sd_proposal)
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
  proposal_value <- rnorm_max_coupling(chain_state1, chain_state2, sd_proposal, sd_proposal)
  proposal1 <- proposal_value[1]
  proposal2 <- proposal_value[2]
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

rinit <- function() rnorm(1, sd = 2)
k <- 5
K <- 20
reject <- TRUE
c_chain <- list()
while(reject){
  c_chain <-  coupled_chains(single_kernel, coupled_kernel, rinit, K = K)
  reject <- (c_chain$meetingtime != 10)
}
c_chain$meetingtime
save(c_chain, file = "diagram.RData")
load("diagram.RData")

tau <- c_chain$meetingtime
tau
chain1 <- c_chain$samples1[,1]
chain2 <- c_chain$samples2[,1]
#
g <- qplot(x = 0:K, y = chain1, geom = "line") + geom_point()
g <- g + geom_line(aes(x = 1:tau, y = chain2[1:tau])) + geom_point(aes(x = 1:tau, y = chain2[1:tau]))
g <- g + xlab("iteration") + ylab("state space") + xlim(-10, K+10)
g <- g + geom_vline(xintercept = tau, linetype = 2)
g <- g + geom_vline(xintercept = k, linetype = 3)
g <- g + geom_vline(xintercept = K, linetype = 3)
g <- g + scale_x_continuous(breaks = c(0, k, tau, K), labels = c("0", TeX("$k = 5$"), TeX("$\\tau = 10$"), TeX("$m = 20$")))
g <- g + geom_segment(data=data.frame(x = (k+1):tau, xend = (k+1):tau, y = chain2[(k+1):tau], yend = chain1[(k+2):(tau+1)]),
                      aes(x = x, xend = xend, y = y, yend = yend), colour = "grey")
g <- g + geom_text(aes(x = -.1, y = -.5, label = paste("X[t]", sep = "")), size = 10, parse = TRUE)
g <- g + geom_text(aes(x = 0.7, y = 1.8, label = paste("Y[t-1]", sep = "")), size = 10, parse = TRUE)
g <- g + geom_text(aes(x = 6, y = chain1[6] - .5*(chain1[6] - chain2[5]), label = paste("Delta[t]", sep = "")), size = 10, parse = TRUE)
g

ggsave(filename = "diagram.pdf", plot = g, width = 15, height = 7)

