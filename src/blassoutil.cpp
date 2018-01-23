#include <RcppEigen.h>
#include "mvnorm.h"
#include "coupling.h"
using namespace Rcpp;
using namespace std;
using namespace Eigen;


//   D_tau_inv <- diag(1/tau2, p, p)
//   A <- XtX + D_tau_inv
//   A_inv <- solve(A)
//   beta <- t(fast_rmvnorm(1, (A_inv %*% XtY)[,1], sigma2 * A_inv))
//   norm <- sum((Y - X %*% beta)^2)
//   betaDbeta <- sum(beta^2 / tau2)
//
// [[Rcpp::export]]
List blassoconditional(const Eigen::VectorXd & Y, const Eigen::MatrixXd & X,
                       const Eigen::VectorXd & XtY, const Eigen::MatrixXd & XtX,
                       const NumericVector tau2, const double sigma2){
  int p = X.cols();

  Eigen::MatrixXd D_tau_inv(p, p);
  D_tau_inv.setIdentity();
  for (int i=0; i < p; i++){
    D_tau_inv(i,i) /= tau2(i);
  }
  Eigen::MatrixXd A_inv(p, p);
  A_inv = (XtX + D_tau_inv).inverse();
  NumericVector mean(wrap(A_inv * XtY));
  NumericMatrix Sigma(wrap(A_inv.array() * sigma2));
  NumericMatrix beta = rmvnorm(1, mean, Sigma);
  Eigen::VectorXd beta_eigen = as<Eigen::VectorXd >(beta);
  double norm = (Y - X * beta_eigen).array().square().sum();
  double betaDbeta = 0;
  for (int i=0; i < p; i++){
    betaDbeta += beta_eigen(i) * beta_eigen(i) / tau2(i);
  }
  return List::create(Named("beta") = beta.row(0), Named("norm") = norm,
                      Named("betaDbeta") = betaDbeta);
}


// [[Rcpp::export]]
List blassoconditional_coupled(const Eigen::VectorXd & Y, const Eigen::MatrixXd & X,
                       const Eigen::VectorXd & XtY, const Eigen::MatrixXd & XtX,
                       const NumericVector & tau21, const NumericVector & tau22, const double sigma21, const double sigma22){
  int p = X.cols();
  Eigen::MatrixXd D_tau_inv1(p, p);
  D_tau_inv1.setIdentity();
  for (int i=0; i < p; i++){
    D_tau_inv1(i,i) /= tau21(i);
  }
  Eigen::MatrixXd A_inv1(p, p);
  A_inv1 = (XtX + D_tau_inv1).inverse();
  NumericVector mean1(wrap(A_inv1 * XtY));
  NumericMatrix Sigma1(wrap(A_inv1.array() * sigma21));
  Eigen::MatrixXd D_tau_inv2(p, p);
  D_tau_inv2.setIdentity();
  for (int i=0; i < p; i++){
    D_tau_inv2(i,i) /= tau22(i);
  }
  Eigen::MatrixXd A_inv2(p, p);
  A_inv2 = (XtX + D_tau_inv2).inverse();
  NumericVector mean2(wrap(A_inv2 * XtY));
  NumericMatrix Sigma2(wrap(A_inv2.array() * sigma22));


  NumericMatrix betas = gaussian_max_couplingC(mean1, mean2, Sigma1, Sigma2);
  // NumericMatrix beta = rmvnorm(1, mean, Sigma);
  Eigen::VectorXd beta_eigen1(p);
  Eigen::VectorXd beta_eigen2(p);
  for (int i=0; i < p; i++){
    beta_eigen1(i) = betas(i,0);
    beta_eigen2(i) = betas(i,1);
  }
  // = as<Eigen::VectorXd >(betas.col(1));
  // = as<Eigen::VectorXd >(betas.col(0));
  double norm1 = (Y - X * beta_eigen1).array().square().sum();
  double norm2 = (Y - X * beta_eigen2).array().square().sum();
  double betaDbeta1 = 0.;
  double betaDbeta2 = 0.;
  for (int i=0; i < p; i++){
    betaDbeta1 += beta_eigen1(i) * beta_eigen1(i) / tau21(i);
    betaDbeta2 += beta_eigen2(i) * beta_eigen2(i) / tau22(i);
  }
  return List::create(Named("beta1") = betas(_,0), Named("beta2") = betas(_,1),
                      Named("norm1") = norm1, Named("norm2") = norm2,
                      Named("betaDbeta1") = betaDbeta1, Named("betaDbeta2") = betaDbeta2);
}
