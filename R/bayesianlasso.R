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
    return(list(chain_state = c(rep(0, p), rep(1, p), 1)))
  }
  #
  gibbs_kernel <- function(state){
    beta <- state$chain_state[1:p]
    tau2 <- state$chain_state[(p+1):(2*p)]
    sigma2 <- state$chain_state[2*p+1]
    res_ <- debiasedmcmc:::blassoconditional(Y, X, XtY, XtX, tau2, sigma2)
    beta <- res_$beta
    norm <- res_$norm
    betaDbeta <- res_$betaDbeta
    sigma2 <- rinversegamma(1, alpha1, 0.5 * (norm + betaDbeta))
    # update tau
    sqrtlambda2sigma2 <- sqrt(lambda2 * sigma2)
    for (component in 1:p){
      tau2[component] <- 1 / rinvgaussian(1, sqrtlambda2sigma2 / abs(beta[component]), lambda2)
    }
    return(list(chain_state = c(beta, tau2, sigma2)))
  }

  coupled_gibbs_kernel <- function(state1, state2){
    beta1 <- state1$chain_state[1:p]
    tau21 <- state1$chain_state[(p+1):(2*p)]
    sigma21 <- state1$chain_state[2*p+1]
    beta2 <- state2$chain_state[1:p]
    tau22 <- state2$chain_state[(p+1):(2*p)]
    sigma22 <- state2$chain_state[2*p+1]
    #
    if (all(tau21 == tau22) && all(sigma21 == sigma22)){
      res_ <- debiasedmcmc:::blassoconditional(Y, X, XtY, XtX, tau21, sigma21)
      beta1 <- res_$beta
      norm1 <- res_$norm
      betaDbeta1 <- res_$betaDbeta
      beta2 <- beta1
      norm2 <- norm1
      betaDbeta2 <- betaDbeta1
    } else {
      res_ <- debiasedmcmc:::blassoconditional_coupled(Y, X, XtY, XtX, tau21, tau22, sigma21, sigma22)
      beta1 <- res_$beta1
      norm1 <- res_$norm1
      betaDbeta1 <- res_$betaDbeta1
      beta2 <- res_$beta2
      norm2 <- res_$norm2
      betaDbeta2 <- res_$betaDbeta2
    }
    sigma2s <- rinversegamma_coupled(alpha1, alpha1,
                               0.5 * (norm1 + betaDbeta1),
                               0.5 * (norm2 + betaDbeta2))
    sigma21 <- sigma2s$xy[1]
    sigma22 <- sigma2s$xy[2]
    # update tau
    sqrtlambda2sigma21 <- sqrt(lambda2 * sigma21)
    sqrtlambda2sigma22 <- sqrt(lambda2 * sigma22)
    for (component in 1:p){
      if (sigma21 == sigma22 && beta1[component] == beta2[component]){
        invtau2 <- rinvgaussian(1, sqrtlambda2sigma21/ abs(beta1[component]), lambda2)
        tau21[component] <- 1 / invtau2
        tau22[component] <- tau21[component]
      } else {
        invtau2s <- rinvgaussian_coupled(sqrtlambda2sigma21 / abs(beta1[component]),
                                         sqrtlambda2sigma22 / abs(beta2[component]),
                                         lambda2, lambda2)
        tau21[component] <- 1 / invtau2s[1]
        tau22[component] <- 1 / invtau2s[2]
      }
    }
    chain_state1 <- c(beta1, tau21, sigma21)
    chain_state2 <- c(beta2, tau22, sigma22)
    identical_ <- all(chain_state1 == chain_state2)
    return(list(state1 = list(chain_state = chain_state1),
                state2 = list(chain_state = chain_state2),
                identical = identical_))
  }
  return(list(rinit = rinit, gibbs_kernel = gibbs_kernel, coupled_gibbs_kernel = coupled_gibbs_kernel))
}
