#include <RcppEigen.h>
#include "mvnorm.h"
using namespace Rcpp;
using namespace std;


// following function returns a n times d matrix
// [[Rcpp::export]]
NumericMatrix fast_rmvnorm_(int nsamples, const NumericVector & mean, const NumericMatrix & covariance){
  int ncols = covariance.cols();
  const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
  Eigen::MatrixXd cholesky_covariance(covariance_.llt().matrixU());
  Eigen::MatrixXd Y(nsamples, ncols);
  for(int j = 0; j < ncols; j++){
    GetRNGstate();
    Y.col(j) = as<Eigen::ArrayXd>(rnorm(nsamples));
    PutRNGstate();
  }
  Y = Y * cholesky_covariance;
  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nsamples; i++){
      Y(i,j) = Y(i,j) + mean(j);
    }
  }
  return wrap(Y);
}

// following function returns a n times d matrix
// [[Rcpp::export]]
NumericMatrix fast_rmvnorm_cholesky_(int nsamples, const NumericVector & mean, const Eigen::MatrixXd & cholesky){
  int ncols = cholesky.cols();
  Eigen::MatrixXd Y(nsamples, ncols);
  // RNGScope scope;
  for(int j = 0; j < ncols; j++){
    GetRNGstate();
    Y.col(j) = as<Eigen::ArrayXd>(rnorm(nsamples));
    PutRNGstate();
  }
  Y = Y * cholesky;
  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nsamples; i++){
      Y(i,j) = Y(i,j) + mean(j);
    }
  }
  return wrap(Y);
}

// [[Rcpp::export]]
NumericVector fast_dmvnorm_(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance){
  const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
  const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::LLT<Eigen::MatrixXd> lltofcov(covariance_);
  Eigen::MatrixXd lower = lltofcov.matrixL();
  Eigen::MatrixXd xcentered(x_);
  double halflogdeterminant = lower.diagonal().array().log().sum();;
  double cst = - (halflogdeterminant) - (x.cols() * 0.9189385);
  for(int j = 0; j < x.cols(); j++){
    for(int i = 0; i < x.rows(); i++){
      xcentered(i,j) = xcentered(i,j) - mean(j);
    }
  }
  Eigen::VectorXd results = -0.5 * lower.triangularView<Eigen::Lower>().solve(xcentered.transpose()).colwise().squaredNorm();
  for (int i = 0; i < results.size(); i++){
    results(i) = results(i) + cst;
  }
  return wrap(results);
}


// [[Rcpp::export]]
NumericVector fast_dmvnorm_cholesky_inverse_(const NumericMatrix & x, const NumericVector & mean, const Eigen::MatrixXd & cholesky_inverse){
  const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::MatrixXd xcentered(x_);
  double halflogdeterminant = cholesky_inverse.diagonal().array().log().sum();;
  double cst = - (-halflogdeterminant) - (x.cols() * 0.9189385);
  for(int j = 0; j < x.cols(); j++){
    for(int i = 0; i < x.rows(); i++){
      xcentered(i,j) = xcentered(i,j) - mean(j);
    }
  }
  Eigen::VectorXd results = -0.5 * (cholesky_inverse.transpose() * xcentered.transpose()).colwise().squaredNorm();
  for (int i = 0; i < results.size(); i++){
    results(i) = results(i) + cst;
  }
  return wrap(results);
}


// returns a list, with 'xy' containing the two generated vectors as columns, and 'identical' (=0 or =1) indicating whether the vectors are equal
// [[Rcpp::export]]
Rcpp::List rmvnorm_max_coupling_(const NumericVector & mu1, const NumericVector & mu2,
                                 const NumericMatrix & Sigma1, const NumericMatrix & Sigma2){
  int p = Sigma1.cols();
  bool identical_ = false;
  NumericMatrix xy(p,2);
  NumericMatrix x1 = fast_rmvnorm_(1, mu1, Sigma1);
  NumericMatrix x2(1,p);
  std::fill(x2.begin(), x2.end(), 0.);
  int nattempts = 0;
  for (int index = 0; index < p; index++){
    xy(index,0) = x1(0,index);
  }
  GetRNGstate();
  NumericVector u = runif(1);
  PutRNGstate();
  double logu = log(u(0));
  if ((fast_dmvnorm_(x1, mu1, Sigma1)(0) + logu) < fast_dmvnorm_(x1, mu2, Sigma2)(0)){
    for (int index = 0; index < p; index++){
      xy(index,1) = xy(index,0);
    }
    identical_ = true;
  } else {
    bool accept = false;
    while (!accept){
      nattempts++;
      x2 = fast_rmvnorm_(1, mu2, Sigma2);
      GetRNGstate();
      u = runif(1);
      PutRNGstate();
      if ((fast_dmvnorm_(x2, mu2, Sigma2)(0) + log(u(0))) > fast_dmvnorm_(x2, mu1, Sigma1)(0)){
        accept = true;
      }
    }
    for (int index = 0; index < p; index++){
      xy(index,1) = x2(0,index);
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical_, Named("nattempts") = nattempts);
}


// [[Rcpp::export]]
Rcpp::List rmvnorm_max_coupling_cholesky(const NumericVector & mu1,
                                         const NumericVector & mu2,
                                         const Eigen::MatrixXd & Cholesky1,
                                         const Eigen::MatrixXd & Cholesky2,
                                         const Eigen::MatrixXd & Cholesky_inverse1,
                                         const Eigen::MatrixXd & Cholesky_inverse2){
  RNGScope scope;
  int p = Cholesky1.cols();
  NumericMatrix xy(p,2);
  NumericMatrix x2(1,p);
  std::fill(x2.begin(), x2.end(), 0.);
  bool identical_ = false;
  NumericMatrix x1 = fast_rmvnorm_cholesky_(1, mu1, Cholesky1);
  for (int index = 0; index < p; index++){
    xy(index,0) = x1(0,index);
  }
  GetRNGstate();
  double u1 = log(runif(1,0,1)(0));
  PutRNGstate();
  if ((fast_dmvnorm_cholesky_inverse_(x1, mu1, Cholesky_inverse1)(0) + u1) <= fast_dmvnorm_cholesky_inverse_(x1, mu2, Cholesky_inverse2)(0)){
    for (int index = 0; index < p; index++){
      xy(index,1) = xy(index,0);
    }
    identical_ = true;
  } else {
    bool accept = false;
    while (!accept){
      x2 = fast_rmvnorm_cholesky_(1, mu2, Cholesky2);
      GetRNGstate();
      double u2 = log(runif(1,0,1)(0));
      PutRNGstate();
      if ((fast_dmvnorm_cholesky_inverse_(x2, mu2, Cholesky_inverse2)(0) + u2) >
            fast_dmvnorm_cholesky_inverse_(x2, mu1, Cholesky_inverse1)(0)){
        accept = true;
      }
    }
    for (int index = 0; index < p; index++){
      xy(index,1) = x2(0,index);
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical_);
}


// Draw from a reflection - maximum coupling of Normal(mu1, Sigma) and Normal(mu2, Sigma)
// As described in "Coupling and Convergence for Hamiltonian Monte Carlo" by Bou-Rabee et al, arXiv:1805.00452v1
// arguments Sigma_chol and inv_Sigma_chol can be obtained e.g. in R as
// chol(Sigma) and solve(chol(Sigma))

// Note that in its current implementation the function does not work for univariate Normals.

// [[Rcpp::export]]
Rcpp::List rmvnorm_reflection_max_coupling_(const Eigen::VectorXd & mu1, const Eigen::VectorXd & mu2,
                                            const Eigen::MatrixXd & Sigma_chol, const Eigen::MatrixXd & inv_Sigma_chol){
  int d = mu1.size();
  NumericMatrix xy(d,2);
  bool identical_ = false;
  // RNGScope scope;
  Eigen::VectorXd scaled_diff = (mu2-mu1).transpose() * inv_Sigma_chol;
  // generate xi and eta, two standard Normal draws
  GetRNGstate();
  Eigen::VectorXd xi = as<Eigen::VectorXd>(rnorm(d, 0., 1.));
  PutRNGstate();
  Eigen::VectorXd eta;
  // define z and e
  Eigen::VectorXd z = - scaled_diff;
  double normz = z.norm();
  if (normz < 1e-15){
    // if normz is very small, most likely the user has supplied identical
    // arguments to mu1 and mu2, which would cause trouble in the forthcoming division by normz
    // thankfully if the user sticks to the function sample_meetingtime, sample_coupled_chains etc
    // then this should never happen.
    identical_ = true;
    xi = mu1 + (xi.transpose() * Sigma_chol).transpose();
    for(int i=0; i<d; i++){
      xy(i,0) = xi(i);
      xy(i,1) = xi(i);
    }
  } else {
    Eigen::ArrayXd e = z.array() / normz;
    GetRNGstate();
    NumericVector utilde = runif(1);
    PutRNGstate();
    double edotxi = (e * xi.array()).sum();
    if (log(utilde(0)) < (-0.5 * (edotxi + normz) * (edotxi + normz) + 0.5 * edotxi * edotxi)){
      eta = (xi.array() + z.array()).matrix();
      identical_ = true;
    } else {
      eta = (xi.array() - 2. * edotxi * e).matrix();
    }
    // construct x ~ Normal(mu1, Sigma) and y ~ Normal(mu2, Sigma) from xi and eta
    xi = mu1 + (xi.transpose() * Sigma_chol).transpose();
    eta = mu2 + (eta.transpose() * Sigma_chol).transpose();
    for(int i=0; i<d; i++){
      xy(i,0) = xi(i);
      if (identical_){
        xy(i,1) = xi(i);
      } else {
        xy(i,1) = eta(i);
      }
    }
  }
  return Rcpp::List::create(Named("xy") = xy, Named("identical") = identical_);
}

