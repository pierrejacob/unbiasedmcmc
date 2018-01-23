# jointly sample ones and zeros
#'@export
coupled_pairs01 <- function(selection1, selection2){
  p <- length(selection1)
  # jointly samples ones
  s1 <- sum(selection1)
  s2 <- sum(selection2)
  w11 <- selection1 / s1
  w12 <- selection2 / s2
  pmin1 <- pmin(w11, w12)
  alpha1 <- sum(pmin1)
  minor1 <- pmin1 / alpha1
  residual1 <- (w11 - pmin1) / (1 - alpha1)
  residual2 <- (w12 - pmin1) / (1 - alpha1)
  ind1 <- c(NA, NA)
  if (runif(1) < alpha1){
    x <- sample(x = 1:p, size = 1, prob = minor1)
    ind1 <- c(x,x)
  } else {
    ind1 <- c(sample(x = 1:p, size = 1, prob = residual1),
              sample(x = 1:p, size = 1, prob = residual2))
  }
  # jointly samples zeros
  w01 <- (1 - selection1) / (p - s1)
  w02 <- (1 - selection2) / (p - s2)
  pmin0 <- pmin(w01, w02)
  alpha0 <- sum(pmin0)
  minor0 <- pmin0 / alpha0
  residual1 <- (w01 - pmin0) / (1 - alpha0)
  residual2 <- (w02 - pmin0) / (1 - alpha0)
  ind0 <- c(NA, NA)
  if (runif(1) < alpha0){
    x <- sample(x = 1:p, size = 1, prob = minor0)
    ind0 <- c(x,x)
  } else {
    ind0 <- c(sample(x = 1:p, size = 1, prob = residual1),
              sample(x = 1:p, size = 1, prob = residual2))
  }
  return(c(ind0, ind1))
}
