#include <Rcpp.h>
#include "RNG.h"
#include "PolyaGamma.h"

using namespace Rcpp;


// [[Rcpp::export]]
double logcosh(double x){
  double result = 0.;
  if (x > 0.){
    result = x + log(1.0 + exp(-2.0*x)) - 0.6931472;
  } else {
    result = -x + log(1.0 + exp(2.0*x)) - 0.6931472;
  }
  return result;
}



// [[Rcpp::export]]
NumericVector xbeta_(const NumericMatrix & X, const NumericVector & beta){
  int n = X.rows();
  NumericVector xbeta(n);
  for (int i = 0; i < n; i++){
    xbeta(i) = 0;
    for (int j = 0; j < X.cols(); j++){
      xbeta(i) = xbeta(i) + X(i,j) * beta(j);
    }
  }
  return(xbeta);
}

double rpg_devroye(int n, double z){
  RNG r;
  PolyaGamma pg(1);
  double x = pg.draw(n, z, r);
  return(x);
}

// [[Rcpp::export]]
NumericMatrix w_rejsamplerC(const NumericVector & beta1,
                            const NumericVector & beta2,
                            const NumericMatrix & X){
  RNGScope scope;
  int n = X.rows();
  NumericMatrix w(n,2);
  NumericVector z1s = abs(xbeta_(X, beta1));
  NumericVector z2s = abs(xbeta_(X, beta2));
  for(int i = 0; i < n; ++i){
    double z1 = z1s(i);
    double z2 = z2s(i);
    double z_min = std::min(z1,z2);
    double z_max = std::max(z1,z2);
    double w_max;
    double w_min = rpg_devroye(1,z_min);
    GetRNGstate();
    double log_u = log(runif(1,0,1)(0));
    PutRNGstate();
    double log_ratio = - 0.5 * w_min * (pow(z_max,2.)-pow(z_min,2.));
    if(log_u < log_ratio){
      w_max = w_min;
    } else {
      w_max = rpg_devroye(1,z_max);
    }
    if(z1<z2){
      w(i,0) = w_min;
      w(i,1) = w_max;
    } else {
      w(i,0) = w_max;
      w(i,1) = w_min;
    }
  }
  return w;
}


// [[Rcpp::export]]
NumericMatrix w_max_couplingC(const NumericVector & beta1,
                              const NumericVector & beta2,
                              const NumericMatrix & X){
  RNGScope scope;
  int n = X.rows();
  NumericMatrix w(n,2);
  NumericVector z1s = abs(xbeta_(X, beta1));
  NumericVector z2s = abs(xbeta_(X, beta2));
  for(int i = 0; i < n; ++i){
    double z1 = z1s(i);
    double z2 = z2s(i);
    double w1 = rpg_devroye(1,z1);
    w(i,0) = w1;
    double w2;
    GetRNGstate();
    double u1 = runif(1)(0);
    PutRNGstate();
    double logaccept1 = logcosh(z2/2.) - 0.5*z2*z2*w1 - (logcosh(z1/2.) - 0.5*z1*z1*w1);
    if(log(u1) <= logaccept1){
      w2 = w1;
    } else {
      bool accept = FALSE;
      while (accept==FALSE){
        w2 = rpg_devroye(1,z2);
        GetRNGstate();
        double u2 = runif(1)(0);
        PutRNGstate();
        double logaccept2 = logcosh(z1/2.) - 0.5*z1*z1*w2 - (logcosh(z2/2.) - 0.5*z2*z2*w2);
        if (log(u2) > logaccept2){
          accept = TRUE;
        }
      }
    }
    w(i,1) = w2;
  }
  return w;
}


