
curve(log(cosh(x)), from = -10, to = 10)


logcosh <- function(x){
  return(abs(x) + log((1 + exp(-2*abs(x)))/2))
}

curve(logcosh, add = TRUE, col = "red", lty = 2)
