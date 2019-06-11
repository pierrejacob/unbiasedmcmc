#' Y and X need to be matrices, and lambda non-negative
#'@export
get_variableselection <- function(Y, X, g, kappa, s0, proportion_singleflip){
  # n <- dim(X)[1]
  p <- dim(X)[2]
  Y2 <- (t(Y) %*% Y)[1,1]
  marginal_likelihood <- function(selection) unbiasedmcmc:::marginal_likelihood_c_2(selection, X, Y, Y2, g)
  prior <- function(selection){
    sumones <- sum(selection)
    ifelse(sumones > s0, -Inf, -kappa * sumones * log(p))
  }
  rinit <- function(){
    x <- rep(0, p)
    x[sample(1:p, min(s0, p), replace = F)] <- (runif(min(s0, p)) < 0.5)
    return(x)
  }
  single_kernel <- function(current_state, current_pdf){
    if (runif(1) < proportion_singleflip){
      # proposal: flip a coordinate
      coordinate <- sample.int(p, 1)
      proposed_state <- current_state
      proposed_state[coordinate] <- 1 - proposed_state[coordinate]
      proposed_prior <- prior(proposed_state)
      if (is.finite(proposed_prior)){
        proposed_pdf <- marginal_likelihood(proposed_state) + proposed_prior
        if (log(runif(1)) < (proposed_pdf - current_pdf)){
          # accept
          current_state <- proposed_state
          current_pdf <- proposed_pdf
        }
      }
    } else {
      # propose to swap a one and a zero
      sumones <- sum(current_state)
      # test if there are ones and zeros to swap
      if (sumones > 0 && sumones < p){
        indices <- unbiasedmcmc:::sample_pair01(current_state)
        proposed_state <- current_state
        proposed_state[indices[1]] <- 1
        proposed_state[indices[2]] <- 0
        proposed_prior <- prior(proposed_state)
        if (is.finite(proposed_prior)){
          proposed_pdf <- marginal_likelihood(proposed_state) + proposed_prior
          if (log(runif(1)) < (proposed_pdf - current_pdf)){
            current_state <- proposed_state
            current_pdf <- proposed_pdf
          }
        }
      } else {
        # if there are no elements to swap
        # stay at the current state
      }
    }
    return(list(state = current_state, pdf = current_pdf))
  }
  # coupled kernel
  coupled_kernel <- function(current_state1, current_state2, current_pdf1, current_pdf2){
    if (runif(1) < proportion_singleflip){
      # propose to flip a coordinate
      coordinate <- sample.int(p, 1)
      proposed_state1 <- current_state1
      proposed_state1[coordinate] <- 1 - proposed_state1[coordinate]
      proposed_prior1 <- prior(proposed_state1)
      proposed_state2 <- current_state2
      proposed_state2[coordinate] <- 1 - proposed_state2[coordinate]
      proposed_prior2 <- prior(proposed_state2)
      logu <- log(runif(1))
      if (is.finite(proposed_prior1)){
        proposed_pdf1 <- marginal_likelihood(proposed_state1) + proposed_prior1
        if (logu < (proposed_pdf1 - current_pdf1)){
          current_state1 <- proposed_state1
          current_pdf1 <- proposed_pdf1
        }
      }
      if (is.finite(proposed_prior2)){
        proposed_pdf2 <- marginal_likelihood(proposed_state2) + proposed_prior2
        if (logu < (proposed_pdf2 - current_pdf2)){
          current_state2 <- proposed_state2
          current_pdf2 <- proposed_pdf2
        }
      }
    } else {
      # propose to swap a one and a zero
      sumones1 <- sum(current_state1)
      sumones2 <- sum(current_state2)
      # test if there are ones and zeros to swap
      if (sumones1 > 0 && sumones1 < p && sumones2 > 0 && sumones2 < p){
        indices <- coupled_pairs01(current_state1, current_state2)
        proposed_state1 <- current_state1
        proposed_state1[indices[1]] <- 1
        proposed_state1[indices[3]] <- 0
        proposed_state2 <- current_state2
        proposed_state2[indices[2]] <- 1
        proposed_state2[indices[4]] <- 0
        proposed_prior1 <- prior(proposed_state1)
        proposed_prior2 <- prior(proposed_state2)
        logu <- log(runif(1))
        if (is.finite(proposed_prior1)){
          proposed_pdf1 <- marginal_likelihood(proposed_state1) + proposed_prior1
          if (logu < (proposed_pdf1 - current_pdf1)){
            current_state1 <- proposed_state1
            current_pdf1 <- proposed_pdf1
          }
        }
        if (is.finite(proposed_prior2)){
          proposed_pdf2 <- marginal_likelihood(proposed_state2) + proposed_prior2
          if (logu < (proposed_pdf2 - current_pdf2)){
            current_state2 <- proposed_state2
            current_pdf2 <- proposed_pdf2
          }
        }
      } else {
        # at least one of the two chain states has no ones or no zeros
        # maybe first chain state has ones and zeros to swap
        if (sumones1 > 0 && sumones1 < p){
          indices <- unbiasedmcmc:::sample_pair01(current_state1)
          proposed_state <- current_state1
          proposed_state[indices[1]] <- 1
          proposed_state[indices[2]] <- 0
          proposed_prior <- prior(proposed_state)
          if (is.finite(proposed_prior)){
            proposed_pdf <- marginal_likelihood(proposed_state) + proposed_prior
            if (log(runif(1)) < (proposed_pdf - current_pdf1)){
              current_state1 <- proposed_state
              current_pdf1 <- proposed_pdf
            }
          }
        }
        # or maybe the second chain has ones and zeros to swap
        if (sumones2 > 0 && sumones2 < p){
          indices <- unbiasedmcmc:::sample_pair01(current_state2)
          proposed_state <- current_state2
          proposed_state[indices[1]] <- 1
          proposed_state[indices[2]] <- 0
          proposed_prior <- prior(proposed_state)
          if (is.finite(proposed_prior)){
            proposed_pdf <- marginal_likelihood(proposed_state) + proposed_prior
            if (log(runif(1)) < (proposed_pdf - current_pdf2)){
              current_state2 <- proposed_state
              current_pdf2 <- proposed_pdf
            }
          }
        }
      }
    }
    return(list(state1 = current_state1, pdf1 = current_pdf1,
                state2 = current_state2, pdf2 = current_pdf2))
  }
  # runs pair of coupled chain for max(tau, m) iterations, and keeping history
  # (thus memory intensive for large p and large m)
  coupled_chains_vs <- function(single_kernel, coupled_kernel, rinit, m = 1, max_iterations = Inf, preallocate = 10){
    chain_state1 <- rinit()
    chain_state2 <- rinit()
    current_pdf1 <- marginal_likelihood(chain_state1) + prior(chain_state1)
    current_pdf2 <- marginal_likelihood(chain_state2) + prior(chain_state2)
    dimstate <- length(chain_state1)
    samples1 <- matrix(nrow = m+preallocate+1, ncol = dimstate)
    nrowsamples1 <- m+preallocate+1
    samples2 <- matrix(nrow = m+preallocate, ncol = dimstate)
    samples1[1,] <- chain_state1
    samples2[1,] <- chain_state2
    current_nsamples1 <- 1
    sres1 <- single_kernel(chain_state1, current_pdf1)
    chain_state1 <- sres1$state
    current_pdf1 <- sres1$pdf
    current_nsamples1 <- current_nsamples1 + 1
    samples1[current_nsamples1,] <- chain_state1
    iter <- 1
    meet <- FALSE
    finished <- FALSE
    meetingtime <- Inf
    while (!finished && iter < max_iterations){
      iter <- iter + 1
      # print(iter)
      if (meet){
        sres1 <- single_kernel(chain_state1, current_pdf1)
        chain_state1 <- sres1$state
        current_pdf1 <- sres1$pdf
        chain_state2 <- chain_state1
      } else {
        # print("coupled step")
        res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
        chain_state1 <- res_coupled_kernel$state1
        current_pdf1 <- res_coupled_kernel$pdf1
        chain_state2 <- res_coupled_kernel$state2
        current_pdf2 <- res_coupled_kernel$pdf2
        if (all(chain_state1 == chain_state2) && !meet){
          # recording meeting time tau
          meet <- TRUE
          meetingtime <- iter
        }
      }
      if ((current_nsamples1+1) > nrowsamples1){
        # print('increase nrow')
        new_rows <- nrowsamples1 - 1
        nrowsamples1 <- nrowsamples1 + new_rows
        samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = dimstate))
        samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = dimstate))
      }
      samples1[current_nsamples1+1,] <- chain_state1
      samples2[current_nsamples1,] <- chain_state2
      current_nsamples1 <- current_nsamples1 + 1
      # stop after max(m, tau) steps
      if (iter >= max(meetingtime, m)){
        finished <- TRUE
      }
    }
    samples1 <- samples1[1:current_nsamples1,,drop=F]
    samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
    return(list(samples1 = samples1, samples2 = samples2,
                meetingtime = meetingtime, iteration = iter, finished = finished))
  }
  # runs pair of coupled chain for max(tau, m) iterations, and only keeping the unbiased
  # estimator of E[h(X)] computed on the fly (thus not very memory intensive)
  unbiasedestimator_vs <- function(single_kernel, coupled_kernel, rinit, h = function(x) x, k = 0, m = 1, max_iterations = Inf){
    chain_state1 <- rinit()
    chain_state2 <- rinit()
    current_pdf1 <- marginal_likelihood(chain_state1) + prior(chain_state1)
    current_pdf2 <- marginal_likelihood(chain_state2) + prior(chain_state2)
    # mcmcestimator computes the sum of h(X_t) for t=k,...,m
    mcmcestimator <- h(chain_state1)
    dimh <- length(mcmcestimator)
    if (k > 0){
      mcmcestimator <- rep(0, dimh)
    }
    # correction computes the sum of min(1, (t - k + 1) / (m - k + 1)) * (h(X_{t+1}) - h(X_t)) for t=k,...,max(m, tau - 1)
    correction <- rep(0, dimh)
    sres1 <- single_kernel(chain_state1, current_pdf1)
    chain_state1 <- sres1$state
    current_pdf1 <- sres1$pdf
    if (k == 0){
      correction <- correction + (min(1, (0 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
    }
    if (k <= 1 && m >= 1){
      mcmcestimator <- mcmcestimator + h(chain_state1)
    }
    iter <- 1
    meet <- FALSE
    finished <- FALSE
    meetingtime <- Inf
    # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
    while (!finished && iter < max_iterations){
      iter <- iter + 1
      if (meet){
        sres1 <- single_kernel(chain_state1, current_pdf1)
        chain_state1 <- sres1$state
        current_pdf1 <- sres1$pdf
        chain_state2 <- chain_state1
        if (k <= iter && iter <= m){
          mcmcestimator <- mcmcestimator + h(chain_state1)
        }
      } else {
        res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
        chain_state1 <- res_coupled_kernel$state1
        current_pdf1 <- res_coupled_kernel$pdf1
        chain_state2 <- res_coupled_kernel$state2
        current_pdf2 <- res_coupled_kernel$pdf2
        if (all(chain_state1 == chain_state2) && !meet){
          # recording meeting time tau
          meet <- TRUE
          meetingtime <- iter
        }
        if (k <= iter){
          if (iter <= m){
            mcmcestimator <- mcmcestimator + h(chain_state1)
          }
          correction <- correction + (min(1, (iter-1 - k + 1)/(m - k + 1))) * (h(chain_state1) - h(chain_state2))
        }
      }
      # stop after max(m, tau) steps
      if (iter >= max(meetingtime, m)){
        finished <- TRUE
      }
    }
    mcmcestimator <- mcmcestimator / (m - k + 1)
    uestimator <- mcmcestimator + correction
    return(list(mcmcestimator = mcmcestimator, correction = correction, uestimator = uestimator,
                meetingtime = meetingtime, iteration = iter, finished = finished))
  }
  return(list(rinit = rinit, marginal_likelihood = marginal_likelihood, prior = prior,
              single_kernel = single_kernel, coupled_kernel = coupled_kernel,
              coupled_chains = coupled_chains_vs, unbiasedestimator = unbiasedestimator_vs))
}
