# load packages
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
# this example is taken from Rosenthal "Parallel Computing and Monte Carlo algorithms", 2000.
# The data are baseball player's batting averages, as in Efron and Morris 1975 or Morris 1983
# In particular, the data are listed in Morris 1983, and available as part of Rosenthal's online code:
# http://probability.ca/jeff/comp/james.c

Y <- c(0.395, 0.375, 0.355, 0.334, 0.313, 0.313, 0.291, 0.269, 0.247, 0.247, 0.224, 0.224,
       0.224, 0.224, 0.224, 0.200, 0.175, 0.148)
ndata <- length(Y)

# The data are Y_i for i in {1,...,K}, with here K = 18.
# The model specifies Y_i ~ Normal(theta_i, V), with V known and fixed to 0.00434
V <- 0.00434
# The prior is theta_i ~ Normal(mu, A), with mu ~ flat prior, and A ~ InverseGamma(a,b)
# with density proportional to exp(-b/x) x^{-a-1}; Rosenthal chooses a = -1, b = 2
a <- -1
b <- 2
# To target the posterior distribution, we follow Rosenthal and consider the following Gibbs sampler.
# A given rest: IG(a + (K-1)/2, b + (sum_i=1^K (theta_i - theta_bar)^2)/2)
# mu given rest: Normal(theta_bar, A / K)
# theta_i given rest: Normal( (mu * V + Y_i * A) / (V + A), A * V / (V + A))
# ... where theta_bar is the average of the theta_i, for i in {1,...,K}


# we store the parameters as (mu, A, theta_1, ..., theta_K) so the parameter space is of dimension 20
# the initialization comes from Rosenthal (except the initialization of mu and A which is irrelevant)
rinit <- function(){
  return(c(0,1,rep(mean(Y), ndata)))
}
# here is the code for the Gibbs sampler
single_kernel <- function(current_state, ...){
  theta <- current_state[3:(ndata+2)]
  theta_bar <- mean(current_state[3:(ndata+2)])
  # update of A given rest
  A <- rigamma(1, a + 0.5 * (ndata-1), b + 0.5 * sum((theta - theta_bar)^2))
  # update of mu given rest
  mu <- rnorm(1, theta_bar, A/ndata)
  # update of each theta_i
  theta <- rnorm(ndata, (mu * V + Y * A) / (V + A), A * V / (V + A))
  return(c(mu, A, theta))
}

niterations <- 5e5
chain <- matrix(nrow = niterations, ncol = ndata+2)
chain[1,] <- rinit()
for (iteration in 2:niterations){
  chain[iteration,] <- single_kernel(chain[iteration-1,])
}
save(niterations, chain, file = "baseball.mcmc.RData")
load(file = "baseball.mcmc.RData")
