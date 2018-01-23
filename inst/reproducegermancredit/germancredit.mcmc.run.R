
# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())

#
##  This example is about the Polya Gamma Gibbs sampler for logistic regression models, as applied to the German credit data of Lichman 2013.

# The data:

# #data <- read.csv('~/Dropbox/PolyaGamma/code/debiasedmcmc/inst/logistic_regression/german_credit.csv')
# data <- read.csv('german_credit.csv')
# Y <- data[,"Creditability"]
# x.categorical <- c('Account.Balance', 'Payment.Status.of.Previous.Credit', 'Purpose', 'Value.Savings.Stocks',
#                    'Length.of.current.employment', 'Sex...Marital.Status', 'Guarantors', 'Most.valuable.available.asset',
#                    'Concurrent.Credits', 'Type.of.apartment', 'Occupation', 'Telephone', 'Foreign.Worker')
# x.quant <- c('Duration.of.Credit..month.', 'Credit.Amount', 'Instalment.per.cent', 'Duration.in.Current.address',
#              'Age..years.', 'No.of.Credits.at.this.Bank', 'No.of.dependents')
# for(x in x.categorical){
#   data[,x] = as.factor(data[,x])
# }
# fmla <- paste('~',paste(c(x.quant,x.categorical),collapse ='+'))
# X <- model.matrix(formula(fmla), data=data)
data(germancredit)
n <- nrow(X)
p <- ncol(X)

# prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
logistic_setting <- logistic_precomputation(Y, X, b, B)

# define MCMC transition kernel
single_kernel <- function(chain_state, logistic_setting){
  zs <- abs(xbeta(logistic_setting$X, t(chain_state)))
  w <- rpg(logistic_setting$n, h=1, z=zs)
  res <- m_and_sigma(w, X, logistic_setting$invB, logistic_setting$KTkappaplusinvBtimesb)
  chain_state <- t(fast_rmvnorm_chol(1, res$m, res$Cholesky))
  return(chain_state)
}

rinit <- function(){
  t(fast_rmvnorm(1, mean = b, covariance = B))
}

niterations <- 10000

chain <- matrix(nrow = niterations, ncol = p)
chain[1,] <- rinit()
for (iteration in 2:niterations){
  chain[iteration,] <- single_kernel(chain[iteration-1,], logistic_setting)
}

save(niterations, chain, file = "germancredit.mcmc.RData")
load("germancredit.mcmc.RData")
hist(chain[,1])
