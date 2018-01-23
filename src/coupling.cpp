#include <Rcpp.h>
#include "mvnorm.h"
#include "RNG.h"
#include "PolyaGamma.h"
#include "coupling.h"


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
NumericMatrix gaussian_max_couplingC(const NumericVector & mu1,
                                     const NumericVector & mu2,
                                     const NumericMatrix & Sigma1,
                                     const NumericMatrix & Sigma2){
    RNGScope scope;
    int p = Sigma1.cols();

    NumericMatrix x(p,2);

    NumericMatrix x1 = rmvnorm(1, mu1, Sigma1);
    x.column(0) = x1(_,0);

    // double d1 = exp(dmvnorm(x1, mu1, Sigma1)(0));
    // double u1 = runif(1,0,d1)(0);
    // double d2 = exp(dmvnorm(x1, mu2, Sigma2)(0));
    double d1 = dmvnorm(x1, mu1, Sigma1)(0);
    double u1 = log(runif(1,0,1)(0));
    double d2 = dmvnorm(x1, mu2, Sigma2)(0);

    // if (u1 <= d2){
    if ((d1 + u1) <= d2){
      x.column(1) = x1(_,0);
    } else {
      bool accept = FALSE;
      while (accept==FALSE){
        NumericMatrix x2 = rmvnorm(1, mu2, Sigma2);
        // double d2 = exp(dmvnorm(x2, mu2, Sigma2)(0));
        // double u2 = runif(1,0,d2)(0);
        // double d1 = exp(dmvnorm(x2, mu1, Sigma1)(0));
        double d2 = dmvnorm(x2, mu2, Sigma2)(0);
        double u2 = log(runif(1,0,1)(0));
        double d1 = dmvnorm(x2, mu1, Sigma1)(0);
        if ((d2 + u2) > d1){
          accept = TRUE;
          x.column(1) = x2(_,0);
        }
      }
    }
    return x;
}

// [[Rcpp::export]]
NumericMatrix gaussian_max_coupling_cholesky(const NumericVector & mu1,
                                     const NumericVector & mu2,
                                     const Eigen::MatrixXd & Cholesky1,
                                     const Eigen::MatrixXd & Cholesky2,
                                     const Eigen::MatrixXd & Cholesky_inverse1,
                                     const Eigen::MatrixXd & Cholesky_inverse2){
  RNGScope scope;
  int p = Cholesky1.cols();

  NumericMatrix x(p,2);

  NumericMatrix x1 = rmvnorm_cholesky(1, mu1, Cholesky1);
  x.column(0) = x1(_,0);

  double d1 = dmvnorm_cholesky_inverse(x1, mu1, Cholesky_inverse1)(0);
  double u1 = log(runif(1,0,1)(0));
  double d2 = dmvnorm_cholesky_inverse(x1, mu2, Cholesky_inverse2)(0);
  if ((d1 + u1) <= d2){
    x.column(1) = x1(_,0);
  } else {
    bool accept = FALSE;
    while (accept==FALSE){
      NumericMatrix x2 = rmvnorm_cholesky(1, mu2, Cholesky2);
      double d2 = dmvnorm_cholesky_inverse(x2, mu2, Cholesky_inverse2)(0);
      double u2 = log(runif(1,0,1)(0));
      double d1 = dmvnorm_cholesky_inverse(x2, mu1, Cholesky_inverse1)(0);
      if ((d2 + u2) > d1){
        accept = TRUE;
        x.column(1) = x2(_,0);
      }
    }
  }
  return x;
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
//  GetRNGstate();
  double x = pg.draw(n, z, r);
//  PutRNGstate();
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

    double log_u = log(runif(1,0,1)(0));
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
    // std::cerr << i << std::endl;
    double z1 = z1s(i);
    double z2 = z2s(i);

    double w1 = rpg_devroye(1,z1);
    w(i,0) = w1;
    // std::cerr << w1 << std::endl;
    // std::cerr << z1 << std::endl;
    double w2;
    double u1 = runif(1)(0);
    double logaccept1 = logcosh(z2/2.) - 0.5*z2*z2*w1 - (logcosh(z1/2.) - 0.5*z1*z1*w1);
    //double upperbound = cosh(z1/2.)*exp(-0.5*z1*z1*w1);
    //std::cerr << upperbound << std::endl;
    //double u1 = runif(1,0,cosh(z1/2.)*exp(-0.5*pow(z1,2.)*w1))(0);
    // if(log(u1) <= cosh(z2/2.)*exp(-0.5*pow(z2,2.)*w1)){
    if(log(u1) <= logaccept1){
      w2 = w1;
    } else {
      // std::cerr << "coupling failed" << std::endl;
      bool accept = FALSE;
      while (accept==FALSE){
        w2 = rpg_devroye(1,z2);
        double u2 = runif(1)(0);
        double logaccept2 = logcosh(z1/2.) - 0.5*z1*z1*w2 - (logcosh(z2/2.) - 0.5*z2*z2*w2);

        // double u2 = runif(1,0,cosh(z2/2.)*exp(-0.5*pow(z2,2.) * w2))(0);
        if (log(u2) > logaccept2){
          accept = TRUE;
        }
      }
    }
    w(i,1) = w2;
  }
  return w;
}


