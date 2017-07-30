# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
setwd("~/Dropbox/PolyaGammaResults/mixture/")

## target distribution
target <- function(x){
  evals <- log(0.5) + dnorm(x, mean = c(-4, 4), sd = 1, log = TRUE)
  return(max(evals) + log(sum(exp(evals - max(evals)))))
}

curve(exp(sapply(x, target)), from = -15, to = 15)

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
Sigma <- diag(sd_proposal^2, 1, 1)
coupled_kernel <- function(chain_state1, chain_state2){
  proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma, Sigma)
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

rinit <- function() rnorm(1, mean = 10)
nsamples <- 10000
# c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#   coupled_chains(single_kernel, coupled_kernel, rinit)
# }
# save(c_chains_, file = "mixture.easy.c_chains.RData")
load(file = "mixture.easy.c_chains.RData")

meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime, nclass = 100)
x <- as.numeric(names(table(meetingtime)))
y <- as.numeric(table(meetingtime)) / nsamples
g <- qplot(x = x, y = 0, yend = y, xend = x, geom = "segment") + xlab("meeting time") + ylab("proportion")
g
# ggsave(filename = "mixture.easy.meetingtime.pdf", plot = g, width = 7, height = 7)

# df <- data.frame()
# for (ichain in 1:2){
#   index <- ((1:nsamples)[meetingtime > 10 & meetingtime < 50])[ichain]
#   c_chain_ <- c_chains_[[index]]
#   df_ <- data.frame(cbind(c_chain_$samples1[2:(c_chain_$iteration+1)], c_chain_$samples2))
#   df_$ichain <- ichain
#   df_$iteration <- 1:c_chain_$iteration
#   df <- rbind(df, df_)
# }
# head(df)
# library(ggthemes)
# ggplot(df, aes(x = iteration, y = X1, colour = factor(ichain), group = ichain)) + geom_line() +
#   geom_line(aes(y = X2)) + scale_color_colorblind()

# index_short <- which(meetingtime>10 & meetingtime <30)[1]
# c_chain_ <- c_chains_[[index_short]]
# matplot(cbind(c_chain_$samples1[2:(c_chain_$iteration+1)], c_chain_$samples2), type = "l", ylab = "x", xlab = "iteration", col = "black")

# index_long <- which(meetingtime>100)[1]
# c_chain_ <- c_chains_[[index_long]]
# matplot(cbind(c_chain_$samples1[2:(c_chain_$iteration+1)], c_chain_$samples2), type = "l", ylab = "x", xlab = "iteration", col = "red")

##
K <- 200
# c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
#   continue_coupled_chains(c_chains_[[irep]], single_kernel, K = K)
# }
# save(c_chains_, c_chains_continued_, file = "mixture.easy.c_chains.RData")
load(file = "mixture.easy.c_chains.RData")

k <- 50
histogram <- histogram_c_chains(c_chains_continued_, 1, k, K, nclass = 100)
g <- plot_histogram(histogram, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "mixture.histogram1.pdf", plot = g, width = 7, height = 7)

# ggsave(filename = "mixture.easy.histogram.pdf", plot = g, width = 7, height = 7)

sd_proposal <- 1
Sigma <- diag(sd_proposal^2, 1, 1)
# c_chains_2 <-  foreach(irep = 1:nsamples) %dorng% {
#   coupled_chains(single_kernel, coupled_kernel, rinit)
# }
# save(c_chains_, c_chains_continued_, c_chains_2, file = "mixture.easy.c_chains.RData")
load(file = "mixture.easy.c_chains.RData")

# c_chains_continued_2 <-  foreach(irep = 1:nsamples) %dorng% {
#   continue_coupled_chains(c_chains_2[[irep]], single_kernel, K = K)
# }
# save(c_chains_, c_chains_continued_, c_chains_2, c_chains_continued_2, file = "mixture.easy.c_chains.RData")
load(file = "mixture.easy.c_chains.RData")

meetingtime2 <- sapply(c_chains_2, function(x) x$meetingtime)
summary(meetingtime2[1:1000])

histogram2 <- histogram_c_chains(c_chains_continued_2[1:1000], 1, k, K, nclass = 100)
g <- plot_histogram(histogram2, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "mixture.histogram2.pdf", plot = g, width = 7, height = 7)

summary(meetingtime2)
histogram2 <- histogram_c_chains(c_chains_continued_2, 1, k, K, nclass = 100)
g <- plot_histogram(histogram2, with_bar = T)
g <- g + stat_function(fun = function(x) sapply(x, function(x_)exp(target(x_))), colour = "red", alpha = 1)
g <- g + xlim(-10,10)
g
ggsave(filename = "mixture.histogram3.pdf", plot = g, width = 7, height = 7)

