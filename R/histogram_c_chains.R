## create histogram
#'@rdname histogram_c_chains
#'@title histogram_c_chains
#'@description Compute histogram approximations of marginal distributions based on coupled Markov chains
#'@param c_chains A list of coupled chains, each as produced by \code{\link{sample_coupled_chains}}
#'@param component An integer specifying which marginal to approximate
#'@param k An integer (see \code{\link{H_bar}})
#'@param m Another integer (see \code{\link{H_bar}})
#'@param breaks A vector indicating how to bin the space (optional)
#'@param nclass An integer specifying the number of bins to aim for, if "breaks" is not specified
#'@param dopar A boolean indicating whether to parallelize the computation (requires doParallel and having registed parallel cores)
#'@export
histogram_c_chains <- function(c_chains, component, k, m, breaks = NULL, nclass = 30, dopar = FALSE){
  nsamples <- length(c_chains)
  lag <- dim(c_chains[[1]]$samples1)[1] - dim(c_chains[[1]]$samples2)[1]
  if (is.null(breaks)){
    breaks <- find_breaks(c_chains, component, nclass, k = k, m = m, lag = lag)
  }
  mids <- create_mids(breaks)
  width <- diff(breaks)[1]
  ### compute histogram
  res_ <- NULL
  if (dopar){
    res_ <- foreach (ibreak = 2:length(breaks), .combine = rbind) %dopar% {
      estimators <- rep(0, nsamples)
      for (irep in 1:nsamples){
        estimators[irep] <- estimator_bin(c_chains[[irep]], component, breaks[ibreak-1], breaks[ibreak], k, m, lag)
      }
      prop <- mean(estimators)
      sd_prop <- sd(estimators) / sqrt(nsamples)
      c(prop, sd_prop)
    }
  } else {
    res_ <- matrix(0, nrow = length(breaks)-1, ncol = 2)
    for (ibreak in 2:length(breaks)){
      estimators <- rep(0, nsamples)
      for (irep in 1:nsamples){
        estimators[irep] <- estimator_bin(c_chains[[irep]], component, breaks[ibreak-1], breaks[ibreak], k, m, lag)
      }
      prop <- mean(estimators)
      sd_prop <- sd(estimators) / sqrt(nsamples)
      res_[ibreak-1,] <- c(prop, sd_prop)
    }
  }
  prop <- res_[,1]
  sd_prop <- res_[,2]
  return(list(mids = mids, breaks = breaks, proportions = prop, sd = sd_prop, width = width))
}

## plot the result of the histogram_c_chains function
#'@export
plot_histogram <- function(histogram, with_bar = TRUE){
  df_ <- data.frame(xmin = histogram$mids - histogram$width/2,
                    xmax = histogram$mids + histogram$width/2,
                    x = histogram$mids,
                    y = histogram$proportions / histogram$width,
                    ymin = (histogram$proportions - 2 * histogram$sd) / histogram$width,
                    ymax = (histogram$proportions + 2 * histogram$sd) / histogram$width)

  g <- ggplot(df_, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + geom_rect(alpha = 0.5)
  if (with_bar){
    g <- g + geom_segment(aes(x = x, xend = x, y = 0, yend = y)) + ylab("density")
  }
  return(g)
}

## create breaks based on list of c_chains
#'@export
find_breaks <- function(c_chains, component, nclass, k, m, lag){
  all_samples <- unlist(sapply(c_chains, function(x) x$samples1[(k+1):(m+1),component]))
  all_samples <- c(all_samples, unlist(sapply(c_chains, function(x) x$samples2[(k):(m+1-lag),component])))
  br <- hist(all_samples, plot=F, nclass = nclass)$breaks
  return(br)
}

#'@export
create_mids <- function(breaks){
  mids <- c()
  for (i in 2:length(breaks)){
    mids <- c(mids, breaks[i-1] + (breaks[i] - breaks[i-1])/2)
  }
  return(mids)
}
## compute estimator of component being in [lower, upper]
#'@export
estimator_bin <- function(c_chains, component, lower, upper, k, m, lag){
  return(estimator_bin_(c_chains, component, lower, upper, k, m, lag))
}
