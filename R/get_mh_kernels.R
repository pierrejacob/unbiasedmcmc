#'@rdname get_mh_kernel
#'@title Get random walk Metropolis-Hastings kernels
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled MH kernels, with Normal random walks.
#' These kernels can then be used in the function \code{\link{coupled_chains}}.
#'@param logtarget function to compute target log-density
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_mh_kernels <- function(logtarget, Sigma_proposal, dimension){
  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)
  # single kernel
  kernel <- function(chain_state, iteration){
    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
    proposal_pdf <- logtarget(proposal_value)
    current_pdf <- logtarget(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(proposal_value)
    } else {
      return(chain_state)
    }
  }
  # coupled kernel
  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)
    current_pdf1 <- logtarget(chain_state1)
    current_pdf2 <- logtarget(chain_state2)
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
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}
