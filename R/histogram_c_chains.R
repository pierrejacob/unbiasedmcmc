## create histogram
#'@export
histogram_c_chains <- function(c_chains, component, k, K, breaks = NULL, nclass = 30){
  nsamples <- length(c_chains)
  if (is.null(breaks)){
    breaks <- find_breaks(c_chains, component, nclass, k)
  }
  mids <- create_mids(breaks)
  width <- diff(breaks)[1]
  # ## compute histogram
  res__ <- foreach (ibreak = 2:length(breaks), .combine = rbind) %dorng% {
    estimators <- rep(0, nsamples)
    for (irep in 1:nsamples){
      estimators[irep] <- estimator_bin_R(c_chains[[irep]], component, breaks[ibreak-1], breaks[ibreak], k, K)
    }
    prop <- mean(estimators)
    sd_prop <- sd(estimators) / sqrt(nsamples)
    c(prop, sd_prop)
  }
  prop <- res__[,1]
  sd_prop <- res__[,2]
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
find_breaks <- function(c_chains, component, nclass, k){
  all_samples <- unlist(sapply(c_chains, function(x) x$samples1[(k+1):nrow(x$samples1),component]))
  all_samples <- c(all_samples, unlist(sapply(c_chains, function(x) x$samples2[k:nrow(x$samples2),component])))
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
estimator_bin_R <- function(c_chains, component, lower, upper, k, K){
  return(estimator_bin(c_chains, component, lower, upper, k, K))
}
