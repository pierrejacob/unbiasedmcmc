#'@rdname get_max_coupling
#'@title Sample from maximally coupled distributions p and q
#'@description Takes two distributions (sampling function and pdf),
#' and return a function to sample from a maximal coupling of these distributions.
#' The function returns a list with
#' "xy": the pair of samples (x,y)
#' "identical": TRUE if x = y, FALSE otherwise
#'@export
get_max_coupling <- function(rp, dp, rq, dq){
  function(){
    x <- rp(1)
    if (dp(x) + log(runif(1)) < dq(x)){
      return(list(xy = c(x,x), identical = TRUE))
    } else {
      reject <- TRUE
      y <- NA
      while (reject){
        y <- rq(1)
        reject <- (dq(y) + log(runif(1)) < dp(y))
      }
      return(list(xy = c(x,y), identical = FALSE))
    }
  }
}
