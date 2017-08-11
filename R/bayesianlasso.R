#' Y and X need to be matrices, and lambda non-negative
#'@export
get_blasso <- function(Y, X, lambda){
  p <- ncol(X)
  n <- nrow(X)
  XtX <- t(X) %*% X
  XtY <- t(X) %*% Y
  alpha1 <- (n-1)/2 + p/2
  lambda2 <- lambda^2
  # naive initialization
  rinit <- function(){
    return(c(rep(0, p), rep(1, p), 1))
  }
  #
  gibbs_kernel <- function(current_state, ...){
    beta <- current_state[1:p]
    tau2 <- current_state[(p+1):(2*p)]
    sigma2 <- current_state[2*p+1]
    D_tau_inv <- diag(1/tau2, p, p)
    A <- XtX + D_tau_inv
    A_inv <- solve(A)
    # update beta
    beta <- t(fast_rmvnorm(1, (A_inv %*% XtY)[,1], sigma2 * A_inv))
    # update sigma
    norm <- sum((Y - X %*% beta)^2)
    betaDbeta <- sum(beta^2 / tau2)
    sigma2 <- rigamma(1, alpha1, 0.5 * (norm + betaDbeta))
    # update tau
    sqrtlambda2sigma2 <- sqrt(lambda2 * sigma2)
    for (component in 1:p){
      tau2[component] <- 1 / rinvgaussian(1, sqrtlambda2sigma2 / abs(beta[component]), lambda2)
    }
    return(c(beta, tau2, sigma2))
  }

  coupled_gibbs_kernel <- function(current_state1, current_state2, ...){
    beta1 <- current_state1[1:p]
    tau21 <- current_state1[(p+1):(2*p)]
    sigma21 <- current_state1[2*p+1]
    beta2 <- current_state2[1:p]
    tau22 <- current_state2[(p+1):(2*p)]
    sigma22 <- current_state2[2*p+1]
    #
    D_tau_inv1 <- diag(1/tau21, p, p)
    D_tau_inv2 <- diag(1/tau22, p, p)
    A1 <- XtX + D_tau_inv1
    A2 <- XtX + D_tau_inv2
    A_inv1 <- solve(A1)
    A_inv2 <- solve(A2)
    # update beta
    if (all(tau21 == tau22) && all(sigma21 == sigma22)){
      mean1 <- A_inv1 %*% XtY
      Sigma1 <- sigma21 * A_inv1
      beta1 <- t(fast_rmvnorm(1, mean1, Sigma1))
      beta2 <- beta1
    } else {
      betas <- gaussian_max_coupling((A_inv1 %*% XtY)[,1], (A_inv2 %*% XtY)[,1], sigma21 * A_inv1, sigma22 * A_inv2)
      beta1 <- betas[,1,drop=F]
      beta2 <- betas[,2,drop=F]
    }
    # update sigma
    norm1 <- sum((Y - X %*% beta1)^2)
    norm2 <- sum((Y - X %*% beta2)^2)
    betaDbeta1 <- sum(beta1^2 / tau21)
    betaDbeta2 <- sum(beta2^2 / tau22)

    sigma2s <- rigamma_coupled(alpha1, alpha1,
                               0.5 * (norm1 + betaDbeta1),
                               0.5 * (norm2 + betaDbeta2))
    sigma21 <- sigma2s[1]
    sigma22 <- sigma2s[2]
    # cat("prop sigma coupled:", mean(sigma21 == sigma22), "\n")
    # update tau
    sqrtlambda2sigma21 <- sqrt(lambda2 * sigma21)
    sqrtlambda2sigma22 <- sqrt(lambda2 * sigma22)
    for (component in 1:p){
      if (sigma21 == sigma22 && beta1[component] == beta2[component]){
        invtau2 <- rinvgaussian(1, sqrtlambda2sigma21/ abs(beta1[component]), lambda2)
        tau21[component] <- 1 / invtau2
        tau22[component] <- tau21[component]
      } else {
        invtau2s <- rinvgaussian_coupled_2(sqrtlambda2sigma21 / abs(beta1[component]),
                                           sqrtlambda2sigma22 / abs(beta2[component]),
                                           lambda2, lambda2)
        tau21[component] <- 1 / invtau2s[1]
        tau22[component] <- 1 / invtau2s[2]
      }
    }
    return(list(chain_state1 = c(beta1, tau21, sigma21), chain_state2 = c(beta2, tau22, sigma22)))
  }
  return(list(rinit = rinit, gibbs_kernel = gibbs_kernel, coupled_gibbs_kernel = coupled_gibbs_kernel))
}
