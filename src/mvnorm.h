#ifndef _INCL_MVNORM_
#define _INCL_MVNORM_
#include <RcppEigen.h>
using namespace Rcpp;

// generate samples from a multivariate Normal distribution
NumericMatrix fast_rmvnorm_(int nsamples, const NumericVector & mean, const NumericMatrix & covariance);
NumericMatrix fast_rmvnorm_cholesky_(int nsamples, const NumericVector & mean, const Eigen::MatrixXd & cholesky);

// evaluate probability density function of a multivariate Normal distribution
NumericVector fast_dmvnorm_(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance);
NumericVector fast_dmvnorm_cholesky_inverse_(const NumericMatrix & x, const NumericVector & mean, const Eigen::MatrixXd & cholesky_inverse);

// Couplings --
// The following functions sample from couplings of multivariate Normals.
// They return a list with two entries:
// "xy" is a matrix with two columns, one for x and one for y, which are Normally distributed vectors.
// "identical" which is a boolean indicating whether x = y.

// sample from maximum coupling of two multivariate Normals
Rcpp::List rmvnorm_max_coupling_(const NumericVector & mu1, const NumericVector & mu2, const NumericMatrix & Sigma1, const NumericMatrix & Sigma2);


 // sample from maximum coupling of two multivariate Normals
Rcpp::List rmvnorm_max_coupling_cholesky(const NumericVector & mu1, const NumericVector & mu2,
                                             const Eigen::MatrixXd & Cholesky1, const Eigen::MatrixXd & Cholesky2,
                                             const Eigen::MatrixXd & Cholesky_inverse1, const Eigen::MatrixXd & Cholesky_inverse2);


// sample from reflection-maximum coupling of two multivariate Normals with common covariance matrix
Rcpp::List rmvnorm_reflection_max_coupling_(const Eigen::VectorXd & mu1, const Eigen::VectorXd & mu2,
                                const Eigen::MatrixXd & Sigma_chol, const Eigen::MatrixXd & inv_Sigma_chol);

#endif

