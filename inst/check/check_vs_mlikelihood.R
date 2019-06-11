# This script attempts at optimizing the code that computes the marginal likelihood
# in the variable selection example, and compares different implementations
library(unbiasedmcmc)
rm(list = ls())
set.seed(1)
#
# simulate data
n <- 100
p <- 5000
SNR <- 3
s_star <- 10
s0 <- 100
sigma0 <- 1
beta_star <- SNR * sqrt(sigma0^2 * log(p) / n) * c(2,-3,2,2,-3,3,-2,3,-2,3, rep(0, p-10))
# independent design
X <- matrix(rnorm(n * p), nrow = n, ncol = p) # fast_rmvnorm_chol(n, rep(0, p), diag(1, p, p))
X <- scale(X)
Y <- X %*% matrix(beta_star, ncol = 1) + rnorm(n, 0, sigma0)
Y <- scale(Y)
Y2 <- (t(Y) %*% Y)[1,1]
g <- p^3
kappa <- 1

rinit <- function(){
  x <- rep(0, p)
  x[sample(1:p, min(s0, p), replace = F)] <- (runif(min(s0, p)) < 0.5)
  return(x)
}

marginal_likelihood <- function(selection){
  Xselected <- as.matrix(X[, selection==1,drop=F])
  if (sum(selection) != 0){
    P <- Xselected %*% tcrossprod(x = solve(crossprod(x = Xselected, y = Xselected)), y = Xselected)
  } else {
    P <- matrix(0, n, n)
  }
  R2_gamma <- t(Y) %*% P %*% Y/ Y2
  return((-sum(selection)/2 * log(1+g) - (n/2) * log(1 + g*(1-R2_gamma)))[1,1])
}

marginal_likelihood_alt <- function(selection){
  Xselected <- as.matrix(X[, selection==1,drop=F])
  if (sum(selection) != 0){
    #P1 <- Xselected %*% solve(t(Xselected) %*% Xselected) %*% t(Xselected)
    # using crossprod, it's a bit faster
    P <- Xselected %*% tcrossprod(x = solve(crossprod(x = Xselected, y = Xselected)), y = Xselected)
  } else {
    P <- matrix(0, n, n)
  }
  posteriorvalue <- -(sum(selection) + 1) / 2 * log(g + 1)
  posteriorvalue <- posteriorvalue - n / 2 * log(Y2 - g / (g + 1) * t(Y) %*% P %*% Y)
  return(posteriorvalue[1,1])
}

cppFunction('
double marginal_likelihood_c_(Eigen::VectorXd selection, const Eigen::MatrixXd & X, const Eigen::VectorXd & Y, double Y2, double g){
  double l = 0.;
  int n = X.rows();
  int p = X.cols();
  int s = selection.sum();
  if (s > 0){
    Eigen::MatrixXd Xselected(n,s);
    int counter = 0;
    for (int column = 0; column < p; column++){
      if (selection(column)){
      Xselected.col(counter) = X.col(column);
      counter++;
      }
    }
    l = Y.transpose() * Xselected * ((Xselected.transpose() * Xselected).inverse()) * Xselected.transpose() * Y;
  } else {
    l = 0.;
  }
  l = -((double) s + 1.) / 2. * log(g + 1.) - (double) n / 2. * log(Y2 - g / (g + 1.) * l);
  return l;
}
            ', depends="RcppEigen")

sel <- rinit()
marginal_likelihood_c_(sel, X, Y, Y2, g)
marginal_likelihood_alt(sel)

marginal_likelihood_c_(rep(0, p), X, Y, Y2, g)
marginal_likelihood_alt(rep(0, p))

selection1 <- rinit()
selection2 <- rinit()

marginal_likelihood(selection1)
marginal_likelihood(selection2)
marginal_likelihood(selection1) - marginal_likelihood(selection2)

marginal_likelihood_alt(selection1)
marginal_likelihood_alt(selection2)
marginal_likelihood_alt(selection1) - marginal_likelihood_alt(selection2)


library(microbenchmark)
microbenchmark(marginal_likelihood(rinit()),
               marginal_likelihood_alt(rinit()),
               marginal_likelihood_c_(rinit(), X, Y, Y2, g))
