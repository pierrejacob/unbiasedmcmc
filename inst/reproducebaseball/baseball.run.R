# load packages
library(unbiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
library(doParallel)
library(doRNG)
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
  return(list(chain_state = c(0,1,rep(mean(Y), ndata))))
}
# here is the code for the Gibbs sampler
single_kernel <- function(state){
  theta <- state$chain_state[3:(ndata+2)]
  theta_bar <- mean(theta)
  # update of A given rest
  A <- rinversegamma(1, a + 0.5 * (ndata-1), b + 0.5 * sum((theta - theta_bar)^2))
  # update of mu given rest
  mu <- rnorm(1, theta_bar, sqrt(A/ndata))
  # update of each theta_i
  theta <- rnorm(ndata, (mu * V + Y * A) / (V + A), sqrt(A * V / (V + A)))
  return(list(chain_state = c(mu, A, theta)))
}

# Now coupled kernel
coupled_kernel <- function(state1, state2){
  theta1 <- state1$chain_state[3:(ndata+2)]
  theta_bar1 <- mean(theta1)
  theta2 <- state2$chain_state[3:(ndata+2)]
  theta_bar2 <- mean(theta2)
  identical_components <- logical(length = ndata+2)
  # update of A given rest
  As <- rinversegamma_coupled(alpha1 = a + 0.5 * (ndata-1), alpha2 = a + 0.5 * (ndata-1),
                              beta1 = b + 0.5 * sum((theta1 - theta_bar1)^2),
                              beta2 = b + 0.5 * sum((theta2 - theta_bar2)^2))
  A1 <- As$xy[1]
  A2 <- As$xy[2]
  identical_components[1] <- As$identical
  # update of mu given rest
  # mus <- gaussian_max_coupling(theta_bar1, theta_bar2, matrix(A1/ndata, 1, 1), matrix(A2/ndata, 1, 1))
  mus <- rnorm_max_coupling(theta_bar1, theta_bar2, sqrt(A1/ndata), sqrt(A2/ndata))
  mu1 <- mus$xy[1]
  mu2 <- mus$xy[2]
  identical_components[2] <- mus$identical
  # # update of each theta_i
  thetas <- matrix(nrow = ndata, ncol = 2)
  for (idata in 1:ndata){
    theta_coupled_ <- rnorm_max_coupling((mu1 * V + Y[idata] * A1) / (V + A1), (mu2 * V + Y[idata] * A2) / (V + A2),
                                         sqrt(A1 * V / (V + A1)), sqrt(A2 * V / (V + A2)))
    thetas[idata,] <- theta_coupled_$xy
    identical_components[2+idata] <- theta_coupled_$identical

  }
  theta1 <- thetas[,1];  theta2 <- thetas[,2]
  return(list(state1 = list(chain_state = c(mu1, A1, theta1)),
              state2 = list(chain_state = c(mu2, A2, theta2)),
              identical = all(identical_components)))
}


### First, get a sample of i.i.d. meeting times
nsamples <- 1000
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_meetingtime(single_kernel, coupled_kernel, rinit)
}
meetingtime.list <- list(c_chains = c_chains_, nsamples = nsamples)
save(meetingtime.list, file = "baseball.tuning.RData")
load(file = "baseball.tuning.RData")
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
table(meetingtime)
# we have seen only times = to 2 and 3
# so we pick k = 4

## Modify nsamples
nsamples <- 1000
m <- 100
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m)
}

save(meetingtime.list, c_chains_, m, file = "baseball.tuning.RData")

# nsamples <- 1000
# m <- 100
# final.c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
#   sample_coupled_chains(single_kernel, coupled_kernel, rinit, m = m)
# }
# save(meetingtime.list, c_chains_, m, final.c_chains_, file = "baseball.tuning.RData")
