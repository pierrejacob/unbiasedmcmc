#'@rdname get_mh_kernel
#'@title Get random walk Metropolis-Hastings kernels
#'@description This function takes a target (specified through its log-pdf)
#' and a covariance matrix for a Normal random walk proposal, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled MH kernels.
#' These kernels can then be used in the functions \code{\link{sample_meetingtime}} or
#' \code{\link{sample_coupled_chains}} or \code{\link{sample_unbiasedestimator}}
#'@param logtarget function taking a vector as input and returning target log-density evaluation
#'@param Sigma_proposal covariance of the Normal random walk proposal
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_mh_kernels <- function(logtarget, Sigma_proposal){
  dimension <- dim(Sigma_proposal)[1]
  Sigma_proposal_chol <- chol(Sigma_proposal)
  Sigma_proposal_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)
  # single kernel
  single_kernel <- function(state){
    chain_state <- state$chain_state
    current_pdf <- state$current_pdf
    proposal_value <- fast_rmvnorm_chol(1, chain_state, Sigma_proposal_chol)
    proposal_pdf <- target(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
  }
  # coupled kernel
  coupled_kernel <- function(state1, state2){
    chain_state1 <- state1$chain_state; current_pdf1 <- state1$current_pdf
    chain_state2 <- state2$chain_state; current_pdf2 <- state2$current_pdf
    proposal_value <- rmvnorm_reflectionmax(chain_state1, chain_state2, Sigma_proposal_chol, Sigma_proposal_chol_inv)
    proposal1 <- proposal_value$xy[,1]; proposal_pdf1 <- target(proposal1)
    if (proposal_value$identical){
      proposal2 <- proposal1; proposal_pdf2 <- proposal_pdf1
    } else {
      proposal2 <- proposal_value$xy[,2]; proposal_pdf2 <- target(proposal2)
    }
    logu <- log(runif(1))
    accept1 <- FALSE; accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
      current_pdf1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      current_pdf2 <- proposal_pdf2
    }
    identical_ <- proposal_value$identical && accept1 && accept2
    return(list(state1 = list(chain_state = chain_state1, current_pdf = current_pdf1),
                state2 = list(chain_state = chain_state2, current_pdf = current_pdf2),
                identical = identical_))
  }
  return(list(single_kernel = single_kernel, coupled_kernel = coupled_kernel))
}
