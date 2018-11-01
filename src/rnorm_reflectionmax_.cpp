#include <RcppEigen.h>
using namespace Rcpp;

// Draw from a reflection - maximum coupling of Normal(mu1, Sigma) and Normal(mu2, Sigma)
// As described in "Coupling and Convergence for Hamiltonian Monte Carlo" by Bou-Rabee et al, arXiv:1805.00452v1
// arguments Sigma_chol and inv_Sigma_chol can be obtained e.g. in R as
// chol(Sigma) and solve(chol(Sigma))
// [[Rcpp::export]]
Rcpp::List rnorm_reflectionmax_(const Eigen::VectorXd & mu1, const Eigen::VectorXd & mu2,
                         const Eigen::MatrixXd & Sigma_chol, const Eigen::MatrixXd & inv_Sigma_chol){
  RNGScope scope;
  int d = mu1.size();
  Eigen::VectorXd scaled_diff = (mu2-mu1).transpose() * inv_Sigma_chol;
  // generate xi and eta, two standard multivariate Normal draws
  Eigen::VectorXd xi = as<Eigen::VectorXd>(rnorm(d, 0., 1.));
  Eigen::VectorXd eta;
  // define z and e
  Eigen::VectorXd z = - scaled_diff;
  double normz = z.norm();
  Eigen::ArrayXd e = z.array() / normz;
  NumericVector utilde = runif(1);
  bool identical = false;
  double edotxi = (e * xi.array()).sum();
  if (log(utilde(0)) < (-0.5 * (edotxi + normz) * (edotxi + normz) + 0.5 * edotxi * edotxi)){
    eta = (xi.array() + z.array()).matrix();
    identical = true;
  } else {
    eta = (xi.array() - 2. * edotxi * e).matrix();
  }
  // construct x ~ Normal(mu1, Sigma) and y ~ Normal(mu2, Sigma) from xi and eta
  xi = mu1 + (xi.transpose() * Sigma_chol).transpose();
  eta = mu2 + (eta.transpose() * Sigma_chol).transpose();
  NumericMatrix xy(d,2);
  for(int i=0; i<d; i++){
    xy(i,0) = xi(i);
    if (identical){
      xy(i,1) = xi(i);
    } else {
      xy(i,1) = eta(i);
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical);
}
