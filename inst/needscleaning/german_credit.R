
#setwd('~/Desktop/temp_unbiased')

# Setup -------------------------------------------------------------------
library(debiasedmcmc)
library(latex2exp)
library(grid)
library(gridExtra)

setmytheme()
rm(list = ls())
set.seed(22)
registerDoParallel(cores = detectCores())

# Read in German credit data ----------------------------------------------

base_dir <- '~/Dropbox'
data <- read.csv(paste0(base_dir,'/PolyaGamma/code/debiasedmcmc/inst/logistic_regression/german_credit.csv'))
Y <- data[,"Creditability"]

x.categorical <- c('Account.Balance', 'Payment.Status.of.Previous.Credit', 'Purpose', 'Value.Savings.Stocks',
                  'Length.of.current.employment', 'Sex...Marital.Status', 'Guarantors', 'Most.valuable.available.asset',
                  'Concurrent.Credits', 'Type.of.apartment', 'Occupation', 'Telephone', 'Foreign.Worker')
x.quant <- c('Duration.of.Credit..month.', 'Credit.Amount', 'Instalment.per.cent', 'Duration.in.Current.address',
             'Age..years.', 'No.of.Credits.at.this.Bank', 'No.of.dependents')
for(x in x.categorical){
  data[,x] = as.factor(data[,x])
}
fmla <- paste('~',paste(c(x.quant,x.categorical),collapse ='+'))
X <- model.matrix(formula(fmla), data=data)

n <- nrow(X)
p <- ncol(X)

# Run GLS logistic regression and set prior ---------------------------------------------------------

fmla2 <- paste0('Y',fmla)
G <- glm(fmla2,data=data)
x = summary(G)$coef
# head(x[order(abs(x[,'t value']),decreasing=TRUE),])
# range(G$coefficients)

# prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
logistic_setting <- logistic_precomputation(Y, X, b, B)

# Set kernels ----------------------------------------------------

# Markov kernel of a single chain
single_kernel <- function(chain_state, logistic_setting){
  zs <- abs(xbeta(logistic_setting$X, t(chain_state)))
  w <- rpg(logistic_setting$n, h=1, z=zs)
  res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(chain_state)
}

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2, logistic_setting, return_ws=FALSE){
  ws <- sample_w(chain_state1, chain_state2, logistic_setting$X)
  betas <- sample_beta(ws$w1, ws$w2, logistic_setting)
  if (all(ws$w1 == ws$w2)){
    betas$beta2 <- betas$beta1
  }
  if(!return_ws){
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2)))
  } else {
    return(list(chain_state1=cbind(betas$beta1), chain_state2=cbind(betas$beta2),w1=ws$w1,w2=ws$w2))
  }
}

# Run PGG and build rinit function --------------------------------------

# PGG estimation
gibbs_iters = 10000
# betas_pg <- pg_gibbs(gibbs_iters, logistic_setting)

#matplot(betas_pg, type = "l", lty = 1)

rinit_prior <- function() t(fast_rmvnorm(1, mean = b, covariance = B))

# Test of the coupled kernel ----------------------------------------------
# niterations <- 100
# current_value1 <- rinit_prior()
# current_value2 <- rinit_prior()
#
# chain1 <- matrix(ncol=p, nrow=niterations)
# chain2 <- matrix(ncol=p, nrow=niterations)
# chain1[1,] <- current_value1
# chain2[1,] <- current_value2
#
# for (t in 2:niterations){
#   current_value <- coupled_kernel(current_value1, current_value2, logistic_setting)
#   chain1[t,] <- current_value1 <- current_value$chain_state1
#   chain2[t,] <- current_value2 <- current_value$chain_state2
# }
#
# distances_ <- sapply(1:niterations, function(t) mean((chain1[t,] - chain2[t,])^2))
# plot(1:niterations, distances_, type = "l")

# Meeting times  ---------------------------------------------
R <- 1000

# c_chains <-  foreach(irep = 1:R) %dorng% {
#   coupled_chains(single_kernel, coupled_kernel, rinit_prior, logistic_setting)
# }
# save(c_chains, file = "pgg.c_chains.RData")
load(file = "pgg.c_chains.RData")

mt_df <- data.frame(meetingtime=sapply(c_chains, function(x) x$meetingtime))
g <- ggplot(data=mt_df,aes(x=meetingtime))+geom_histogram(aes(y = ..density..), binwidth=10)
g <- g + xlab("meeting time")
g
ggsave(filename="pgg.meetingtime.pdf", plot=g, width=7, height=7)

summary(mt_df)

# Test of the coupled kernel ----------------------------------------------
#
niterations <- 1000
current_value1 <- rinit_prior()
current_value2 <- rinit_prior()

chain1_w <- chain2_w <- matrix(ncol=n, nrow=niterations)
chain1_beta <- chain2_beta <- matrix(ncol=p, nrow=niterations)
chain1_beta[1,] <- current_value1
chain2_beta[1,] <- current_value2

for(t in 2:1000){
  current_value <- coupled_kernel(current_value1, current_value2, logistic_setting, return_ws=TRUE)
  chain1_w[t-1,] <- current_value$w1
  chain2_w[t-1,] <- current_value$w2
  chain1_beta[t,] <- current_value1 <- current_value$chain_state1
  chain2_beta[t,] <- current_value2 <- current_value$chain_state2
  if(all(current_value1==current_value2)) break
}
meetingtime=t
# t

#trace of w1 and w2
for(idx in 1:1000){
  w1 <- chain1_w[1:(meetingtime-1),idx]
  w2 <- chain2_w[1:(meetingtime-1),idx]
  if(mean(w1==w2)<.6) break
}
idx
w12df <- data.frame(iter=1:(meetingtime-1),chain='X',value=w1)
w12df <- rbind(w12df,data.frame(iter=1:(meetingtime-1),chain='Y',value=w2))

g0 <- ggplot(w12df,aes(x=iter,y=value, col=chain)) + geom_line() + xlab('iteration')  + ylab(TeX('trace of w_{276}'))
g0 <- g0 + scale_colour_grey(start=.75,end=0) + theme(legend.position=c(0.1, .9))
g0

# compute distances and #met for plotting
nmet_w <- sapply(1:meetingtime, function(t) sum(chain1_w[t,]==chain2_w[t,]))
nmet_w[meetingtime] <- n
distances_w <- sapply(1:meetingtime, function(t) sqrt(sum((chain1_w[t,] - chain2_w[t,])^2)))
distances_w[meetingtime] <- 0
distances_beta <- sapply(1:meetingtime, function(t) sqrt(sum((chain1_beta[t,] - chain2_beta[t,])^2)))
df_beta <- data.frame(iter=1:meetingtime, dist_beta = distances_beta[1:meetingtime])
df_w <- data.frame(iter=1:meetingtime, dist_w = distances_w[1:meetingtime], nmet_w = nmet_w[1:meetingtime])

nmet_w <- sapply(1:niterations, function(t) sum(chain1_w[t,]==chain2_w[t,]))
nmet_w[meetingtime] <- n
distances_w <- sapply(1:niterations, function(t) sqrt(sum((chain1_w[t,] - chain2_w[t,])^2)))
distances_w[meetingtime] <- 0
distances_beta <- sapply(1:niterations, function(t) sqrt(sum((chain1_beta[t,] - chain2_beta[t,])^2)))
df_beta <- data.frame(iter=1:meetingtime, dist_beta = distances_beta[1:meetingtime])
df_w <- data.frame(iter=1:meetingtime, dist_w = distances_w[1:meetingtime], nmet_w = nmet_w[1:meetingtime])

g1 <- ggplot(df_beta,aes(x=iter,y=dist_beta))+geom_line() + xlab('iteration') + ylab('distance')
  #ylab(TeX("$||\\beta_1-\\beta_2||_2$"))
g1
ggsave(filename = "pgg.distbeta.pdf", plot = g1, width = 7, height = 7)

g2 <- ggplot(df_w,aes(x=iter,y=dist_w))+geom_line() + xlab('iteration') + ylab('distance')
#+ ylab(TeX("$||w_1-w_2||_2$"))
g2
ggsave(filename = "pgg.distw.pdf", plot = g2, width = 7, height = 7)

g3 <- ggplot(df_w,aes(x=iter,y=nmet_w))+geom_line() + xlab('iteration') + ylab("count")
g3
ggsave(filename = "pgg.nmetw.pdf", plot = g3, width = 7, height = 7)

#grid.arrange(g0, g3, g2, g1, ncol = 2, nrow=2)



# Extend runs -------------------------------------------------------------

max_t = max(sapply(c_chains, function(x) x$iteration))

c_chains_continued <- foreach(chain = c_chains) %dorng% {
  continue_coupled_chains(chain, single_kernel, 2*max_t, logistic_setting)
}

# c_chains_continued <- foreach(chain = c_chains_continued) %dorng% {
#   continue_coupled_chains(chain, single_kernel, 2*max_t, logistic_setting)
# }

max_t = max(sapply(c_chains_continued, function(x) x$iteration))
observed_tau = sapply(c_chains_continued, function(x) x$meetingtime)

save(c_chains, c_chains_continued, file = "pgg.c_chains.RData")
load(file = "pgg.c_chains.RData")


# Find optimal k ----------------------------------------------------------

idx1 <- which('Instalment.per.cent' == colnames(X))
idx2 <- which('Duration.in.Current.address' == colnames(X))
idx <- idx2
colnames(X)[idx]

#k_vec2 <- seq(0,max_t,by=10)
k_vec2 <- seq(0,250,by=5)
lkv = length(k_vec2)
ae_vec2 = rep(NA,lkv)
for(i in 1:lkv){
  k = k_vec2[[i]]
  print(k)
  mean_estimators1 <- foreach(chain = c_chains_continued) %do% {
    H_bar(chain, h = function(x) x[[idx]], k=k, K=k)
  }
  v <- var(as.numeric(mean_estimators1))
  e_max <- mean(pmax(observed_tau,k))
  ae_vec2[[i]] = 1/(v*e_max)
}

df2 <- data.frame(k=k_vec2,ae=ae_vec2)
df2 <- df2[df2$k <= 250,]
k_opt <- k_vec2[which.max(ae_vec2)]
g<- ggplot(df2,aes(x=k,y=ae)) + geom_line() + geom_point() + ylab('asymptotic efficiency')
g

ggsave(filename = "pgg.asympt_eff.pdf", plot = g, width = 7, height = 7)

k_opt
mean(observed_tau<k_opt)

K_opt <- k_opt * 10

k_opt
K_opt

# Create histograms for two variable estimates ---------------------------------

c_chains_continued2 <- foreach(chain = c_chains_continued) %dorng% {
  continue_coupled_chains(chain, single_kernel, K_opt, logistic_setting)
}

save(c_chains, c_chains_continued, c_chains_continued2, file = "pgg.c_chains.RData")
load(file = "pgg.c_chains.RData")


# histogram
idx1 <- which('Instalment.per.cent' == colnames(X))
idx2 <- which('Duration.in.Current.address' == colnames(X))
idx <- idx1

nclass <- 27 #27 for idx1, 24 for idx2
rng <- range(find_breaks(c_chains_continued2, idx, nclass, k_opt))
breaks = seq(rng[1],rng[2],length=nclass)

# histogram1 <- histogram_c_chains(c_chains_continued2, idx, k_opt, K_opt, nclass = nclass)
histogram1 <- histogram_c_chains(c_chains_continued2, idx, k_opt, K_opt, breaks=breaks)

# niterations <- nrow(betas_pg)
# hist_mcmc <- hist(betas_pg[1000:niterations,idx], breaks = histogram1$breaks, plot = FALSE)
#
# g1 <- plot_histogram(histogram1, with_bar = TRUE) + xlab(TeX("$\\beta$")) + ylab("density")
# g1 <- g1 + geom_line(aes(x = hist_mcmc$mids, y = hist_mcmc$density), colour = "red")
# g1
#
# ggsave(filename = "pgg.histogram1.pdf", plot = g1, width = 7, height = 7)






