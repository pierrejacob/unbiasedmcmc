#include <Rcpp.h>
#include "mvnorm.h"
#include "RNG.h"
#include "PolyaGamma.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gaussian_max_couplingC(const NumericVector & mu1,
                                     const NumericVector & mu2,
                                     const NumericMatrix & Sigma1,
                                     const NumericMatrix & Sigma2){
    RNGScope scope;
    int p = Sigma1.cols();

    NumericMatrix x(p,2);

    NumericMatrix x1 = rmvnorm(1, mu1, Sigma1);
    x.column(0) = x1(_,0);

    double d1 = exp(dmvnorm(x1, mu1, Sigma1)(0));
    double u1 = runif(1,0,d1)(0);
    double d2 = exp(dmvnorm(x1, mu2, Sigma2)(0));
    if (u1 <= d2){
      x.column(1) = x1(_,0);
    } else {
      bool accept = FALSE;
      while (accept==FALSE){
        NumericMatrix x2 = rmvnorm(1, mu2, Sigma2);
        double d2 = exp(dmvnorm(x2, mu2, Sigma2)(0));
        double u2 = runif(1,0,d2)(0);
        double d1 = exp(dmvnorm(x2, mu1, Sigma1)(0));
        if (u2 > d1){
          accept = TRUE;
          x.column(1) = x2(_,0);
        }
      }
    }
    return x;
}

// [[Rcpp::export]]
NumericVector gaussian_max_coupling_cholesky(const NumericVector & mu1,
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

  double d1 = exp(dmvnorm_cholesky_inverse(x1, mu1, Cholesky_inverse1)(0));
  double u1 = runif(1,0,d1)(0);
  double d2 = exp(dmvnorm_cholesky_inverse(x1, mu2, Cholesky_inverse2)(0));
  if (u1 <= d2){
    x.column(1) = x1(_,0);
  } else {
    bool accept = FALSE;
    while (accept==FALSE){
      NumericMatrix x2 = rmvnorm_cholesky(1, mu2, Cholesky2);
      double d2 = exp(dmvnorm_cholesky_inverse(x2, mu2, Cholesky_inverse2)(0));
      double u2 = runif(1,0,d2)(0);
      double d1 = exp(dmvnorm_cholesky_inverse(x2, mu1, Cholesky_inverse1)(0));
      if (u2 > d1){
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

//
// w_rejsampler <- function(beta1, beta2, X, max_iter=1000){
// # we can do rejection sampling, sampling from z1 and aiming for z2,
// # provided z2 > z1. The ratio of densities is proportional to
//   logratio <- function(omega, z_min, z_max){
//     return(-omega * 0.5 * (z_max^2 - z_min^2))
//   }
//
//   w1 <- rep(0., n)
//   w2 <- rep(0., n)
//   z1s <- abs(xbeta(X, beta1))
//   z2s <- abs(xbeta(X, beta2))
//
//   for (i in 1:n){
//     z_i <- c(z1s[i], z2s[i])
//     z_min <- min(z_i)
//     z_max <- max(z_i)
//
//     proposal <- rpg(num = 1, h = 1, z = z_min)
//     w_min <- proposal
//
//     iter <- 0
//     if (log(runif(1)) < logratio(proposal, z_min, z_max)){
//       w_max <- proposal
//     } else {
//       w_max <- rpg(num = 1, h = 1, z = z_max)
//     }
//
//     if (which.min(z_i) == 1){
//       w1[i] <- w_min
//       w2[i] <- w_max
//     } else {
//       w2[i] <- w_min
//       w1[i] <- w_max
//     }
//   }
//   return(list(w1 = w1, w2 = w2))
// }
//

