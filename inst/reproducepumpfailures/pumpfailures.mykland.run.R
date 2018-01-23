# load packages

##  This example is about failures of nuclear pumps. It's classic (e.g. Example 10.17 in Robert & Casella Monte Carlo Statistical Methods)
## It's used as an example in Murdoch and Green's perfect samplers paper and also Reutter and Johnson 1995
## about using coupled chains to monitor MCMC convergence

# The data:
# number of failures
s <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
# times
t <- c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
# the model says s_k ~ Poisson(lambda_k * t_k), for k = 1,...,10
ndata <- 10
# and lambda_k ~ Gamma(alpha,beta), beta ~ Gamma(gamma, delta)
alpha <- 1.802
gamma <- 0.01
delta <- 1
# full conditionasl:
# lambda_k given rest: Gamma(alpha + s_k, beta + t_k)
# beta given rest: Gamma(gamma + 10*alpha, delta + sum_{k=1}^10 lambda_k)

p_D <- function(Lambda,d1,d2){
  F_d1 <- pgamma(d1,gamma+10*alpha,Lambda+delta)
  F_d2 <- pgamma(d2,gamma+10*alpha,Lambda+delta)
  return(F_d2-F_d1)
}

d_function <- function(Lambda, Lambda_tilde, d1, d2){
    if(Lambda<Lambda_tilde){
      return(d1)
    } else {
      return(d2)
    }
}

s_function <- function(Lambda, Lambda_tilde, d1, d2){
  s <- p_D(Lambda,d1,d2)
  s <- s * ((Lambda+delta)/(Lambda_tilde+delta))^(gamma+10*alpha)
  s <- s * exp((Lambda_tilde-Lambda)*d_function(Lambda, Lambda_tilde, d1, d2))
  return(s)
}

r_function <- function(x, y, Lambda_tilde=6.7, d1=2.35-1.1*.69, d2=2.35+1.1*.69){
  lambda_x <- x[1:ndata]
  Lambda_x = sum(lambda_x)
  beta_y <-   x[ndata+1]

  if(d1 <= beta_y & beta_y <= d2){
    r <- exp((Lambda_tilde - Lambda_x)*
               (d_function(Lambda_x,Lambda_tilde,d1,d2) - beta_y))
    return(r)
  } else {
    return(0)
  }
}

# run split Gibbs sampler

single_kernel <- function(current_state){
  lambda <- current_state[1:ndata]
  beta <- current_state[ndata+1]
  for (k in 1:ndata){
    lambda[k] <- rgamma(1, shape = alpha + s[k], rate = beta + t[k])
  }
  beta <- rgamma(1, shape = gamma + 10*alpha, rate = delta + sum(lambda))
  return(c(lambda, beta))
}

rinit <- function(){
  return(rep(1, ndata+1))
}

n_rep = 1000

nt_vec = rep(NA,n_rep)
tl_vec = c()
beta_hat_vec = rep(NA,n_rep)

for(jj in 1:n_rep){
  print(jj)
  niterations <- 5000
  chain <- matrix(nrow = niterations, ncol = ndata+1)
  chain[1,] <- rinit()
  for (iteration in 2:niterations){
    chain[iteration,] <- single_kernel(chain[iteration-1,])
  }

  S_vec = rep(0, niterations)
  for(iteration in 1:(niterations-1)){
    r_xy = r_function(chain[iteration,],chain[iteration+1,])
    S_vec[[iteration]] <- rbinom(1, 1, r_xy)
  }

  T_vec <- which(S_vec==1)
  n_tours <- length(T_vec)-1
  t_length = T_vec[2:length(T_vec)] - T_vec[1:(length(T_vec)-1)]

  tl_vec = c(tl_vec,t_length)
  nt_vec[jj] = n_tours

  wh0 = T_vec[[1]]+1
  wh1 = T_vec[[length(T_vec)]]
  est_n = sum(chain[wh0:wh1,ncol(chain),drop=FALSE])
  est_d = wh1-wh0+1
  beta_hat_vec[jj] = est_n/est_d
}

mean(tl_vec)
mean(nt_vec)
mean(beta_hat_vec)
format(var(beta_hat_vec),scientific = TRUE)
(5000*var(beta_hat_vec))^(-1)


Y_chain <- matrix(nrow = n_tours, ncol = ndata+1)
N_vec <- rep(NA, n_tours)

#SRQ plot
#plot((1:(n_tours+1))/(n_tours+1),T_vec/tail(T_vec,1),type='l')

# estimators
wh_cv = c()
for(i in 1:n_tours){
  wh <- (T_vec[[i]]+1):T_vec[[i+1]]
  wh_cv = c(wh_cv,wh)
  Y_chain[i,] <- colSums(chain[wh,,drop=FALSE])
  N_vec[[i]] <- length(wh)
}

#length(wh0:wh1)
#length(wh_cv)

theta_hat <- colSums(Y_chain)/sum(N_vec)
sigma2_hat <- rep(NA, ndata+1)
for(i in 1:(ndata+1)){
  sigma2_hat[[i]] = mean((Y_chain[,i] - theta_hat[i]*N_vec)^2)/mean(N_vec)^2
}
cbind(theta_hat,sigma2_hat)

