#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rinvgaussian_c(int n, double mu, double lambda){
  RNGScope scope;
  NumericVector results(n);
  NumericVector nu = rnorm(n);
  NumericVector z = runif(n);
  NumericVector x = nu * nu;
  x = mu + mu*mu * x / (2. * lambda) - mu / (2. * lambda) * sqrt(4. * mu * lambda * x + mu*mu * x*x);
  for (int i = 0; i < n; i++){
    if (z(i) <= mu / (mu + x(i))){
      results(i) = x(i);
    } else {
      results(i) = mu * mu / x(i);
    }
  }
  return results;
}

double dinvgaussian_c(double x, double mu, double lambda){
  return 0.5 * log(lambda/6.283185) - 1.5 * log(x) - lambda * (x - mu) * (x - mu) / (2 * mu * mu * x);
}


// [[Rcpp::export]]
NumericVector rinvgaussian_coupled_c(double mu1, double mu2, double lambda1, double lambda2){
  RNGScope scope;
  NumericVector results(2);
  NumericVector x = rinvgaussian_c(1, mu1, lambda1);
  NumericVector u = runif(1);
    if (x(0) < 1e-20){
      x(0) = 1e-20;
    }
    if (dinvgaussian_c(x(0), mu1, lambda1) + log(u(0)) < dinvgaussian_c(x(0), mu2, lambda2)){
      results(0) = x(0);
      results(1) = x(0);
    } else {
      bool reject = true;
      NumericVector y(1);
      while (reject){
        y = rinvgaussian_c(1, mu2, lambda2);
        if (y(0) < 1e-20){
          y(0) = 1e-20;
        }
        u = runif(1);
        reject = (dinvgaussian_c(y(0), mu2, lambda2) + log(u(0)) < dinvgaussian_c(y(0), mu1, lambda1));
      }
      results(0) = x(0);
      results(1) = y(0);
    }

  return results;
}

