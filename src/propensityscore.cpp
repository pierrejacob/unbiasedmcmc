#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector beta2e_(const NumericVector & beta, const NumericMatrix & C){
  double betatimescovar;
  int nbeta = beta.size();
  int nobservations = C.nrow();
  NumericVector e(nobservations);
  for (int i = 0; i < nobservations; i++){
    betatimescovar = beta(0);
    for (int j = 0; j < nbeta-1; j++){
      betatimescovar += beta(j+1) * C(i,j);
    }
    e(i) = 1. / (1. + exp(-betatimescovar));
  }
  return e;
}

// function that takes a vector of doubles
// compute the quintiles
// and return a vector of indicators of quintile memberships,
// ie with values in {1,2,3,4,5}.
// [[Rcpp::export]]
IntegerVector cut_in_fifth_(const NumericVector & x){
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  NumericVector q = NumericVector::create(0.2, 0.4, 0.6, 0.8);
  int s = x.size();
  int nqs = q.size();
  NumericVector quintiles(nqs);
  for (int i = 0; i < nqs; i ++){
    quintiles(i) = y[(int) (s * q[i])];
  }
  IntegerVector indicators(s);
  int indicator;
  for (int i = 0; i < s; i ++){
    indicator = 0;
    while ((x(i) >= quintiles(indicator)) && (indicator < (nqs - 1))){
      indicator ++;
    }
    if (x(i) >= quintiles(indicator)){ indicator ++; }
    indicators(i) = indicator + 1; // so that the result is between 1 and 5.
  }
  return indicators;
}


// [[Rcpp::export]]
NumericVector propensity_module2_loglik2_(NumericMatrix theta1s, NumericMatrix theta2s, const NumericVector & X, const NumericMatrix & C, const NumericVector & Y){
  int n = theta1s.nrow();
  NumericVector evals(n);
  std::fill(evals.begin(), evals.end(), 0);
  int theta_dim = 7;
  double logpr, betatimescovar;
  for (int itheta = 0; itheta < n; itheta++){
    NumericVector beta = theta1s.row(itheta);
    int nobservations = X.size();
    NumericVector e = beta2e_(beta, C);
    IntegerVector indicators = cut_in_fifth_(e);
    NumericVector theta = theta2s.row(itheta);
    double z;
    for (int i = 0; i < nobservations; i++){
      z = theta(0) + theta(1) * X(i);
      if (indicators(i) > 1){
        z = z + theta(0 + indicators(i));
      }
      logpr = - log(1 + exp(-z));
      evals(itheta) += logpr + (1 - Y(i)) * (-z);
    }
  }
  return evals;
}
