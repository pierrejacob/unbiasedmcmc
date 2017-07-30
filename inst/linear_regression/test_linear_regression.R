

# Setup -------------------------------------------------------------------
library(debiasedmcmc)
rm(list = ls())
setmytheme()
set.seed(21)
registerDoParallel(cores = detectCores())


# Draw input data ---------------------------------------------------------
n <- 1000
p <- 25

# covariates
meanX <- rep(0, p)
sigmaX <- diag(1, nrow = p, ncol = p)
X <- fast_rmvnorm(n, meanX, sigmaX)

# outcome
sparsity <- ceiling(p/4)
theta_star <- c(rep(1,sparsity),rep(0,p-sparsity))
theta_star <- cbind(theta_star)
sigma_Y <- 2.5
Y <- X %*% theta_star + rnorm(n,sd=sigma_Y)

# prior
b <- matrix(0, nrow = p, ncol = 1)
tau <- 10

# not clear yet if we'll need this
linear_precomputation <- function(Y, X, sigma_Y, b, tau){
  out <- vector(p,mode='list')
  for(j in 1:p){
    Xj <- X[,j,drop=FALSE]
    Xnotj <- X[,-j,drop=FALSE]
    theta <- theta_star

    theta_notj <- theta[-j]

    XtXinv <- solve(crossprod(Xj))
    regrcoef_Y_on_Xj <- XtXinv %*% crossprod(Xj,Y)
    regrcoefs_Xnotj_on_Xj <- XtXinv %*% crossprod(Xj,Xnotj)

    Sigma_data <- drop(XtXinv) * sigma_Y
    Sigma_gibbs <- 1/(1/Sigma_data + 1/tau) # output this for each j
    wt_data <- drop((1/Sigma_data)/(1/Sigma_data + 1/tau))

    out[[j]] <- list(j=j,Sigma_gibbs=Sigma_gibbs,wt_data=wt_data,
                     regrcoef_Y_on_Xj=regrcoef_Y_on_Xj,
                     regrcoefs_Xnotj_on_Xj=regrcoefs_Xnotj_on_Xj)
  }
  return(list(n=nrow(X), p=ncol(X), X=X, Y=Y, sigma_Y=sigma_Y, tau=tau, b=b, matrices=out))
}

linear_setting <- linear_precomputation(Y, X, sigma_Y, b, tau)

conditional_moments <- function(j, theta, linear_setting){
  mtcs <- linear_setting$matrices[[j]]
  Sigma <- matrix(mtcs$Sigma_gibbs,1,1)
  mu_data <- mtcs$regrcoef_Y_on_Xj - mtcs$regrcoefs_Xnotj_on_Xj %*% theta[-j,1,drop=FALSE]
  mu <- drop(wt_data * mu_data + (1-wt_data) * linear_setting$b[j])
  return(list(mu=mu,Sigma=Sigma))
}


# Define Markov kernels for 1x and 2x chains ------------------------------

rinit_prior <- function() t(fast_rmvnorm(1, mean = b, covariance = diag(tau,p)))

# Markov kernel of a single chain
single_kernel <- function(chain_state, linear_setting){
  for(j in 1:linear_setting$p){
    cm <- conditional_moments(j, chain_state, linear_setting)
    chain_state[[j]] <- rnorm(1, cm$mu, cm$Sigma)
  }
  return(chain_state)
}

# Markov kernel of the coupled chain
coupled_kernel <- function(chain_state1, chain_state2, linear_setting){
  for(j in 1:linear_setting$p){
    cm1 <- conditional_moments(j, chain_state1, linear_setting)
    cm2 <- conditional_moments(j, chain_state2, linear_setting)

    gmc <- gaussian_max_coupling(cm1$mu, cm2$mu, cm1$Sigma, cm2$Sigma)
    chain_state1[j,1] <- gmc[1]
    chain_state2[j,1] <- gmc[2]
  }

  return(list(chain_state1=chain_state1,
              chain_state2=chain_state2))
}


### test of the coupled kernel
niterations <- 1000
current_value1 <- rinit_prior()
current_value2 <- rinit_prior()
chain1 <- matrix(ncol=p, nrow=niterations)
chain2 <- matrix(ncol=p, nrow=niterations)
chain1[1,] <- current_value1
chain2[1,] <- current_value2
for (t in 2:niterations){
  current_value <- coupled_kernel(current_value1, current_value2, linear_setting)
  current_value1 <- current_value$chain_state1
  current_value2 <- current_value$chain_state2
  chain1[t,] <- current_value1
  chain2[t,] <- current_value2
}

distances_ <- sapply(1:niterations, function(t) mean((chain1[t,] - chain2[t,])^2))
mean(diff(chain1[,1]) > 1e-10)
plot(1:niterations, distances_, type = "l")

it <- floor(seq(from = 1, to  = niterations, length.out = 1000))
qplot(x=it, y=chain1[it,1], geom = "line") + ylab("X") + xlab("iteration") + geom_line(aes(y = chain2[it,1]), colour = "red")
qplot(x=it, y=chain1[it,2], geom = "line") + ylab("X") + xlab("iteration") + geom_line(aes(y = chain2[it,2]), colour = "red")
acf(chain1[,1])

#
# NOTE: we're getting an error here. Need to debug this
nsamples <- 10
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(single_kernel, coupled_kernel, rinit_prior, linear_setting)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)



