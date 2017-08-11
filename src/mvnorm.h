#ifndef _INCL_MVNORM_
#define _INCL_MVNORM_
#include <RcppEigen.h>
using namespace Rcpp;

// generate samples from a multivariate normal distribution
NumericMatrix rmvnorm(int nsamples, const NumericVector & mean, const NumericMatrix & covariance);
NumericMatrix rmvnorm_cholesky(int nsamples, const NumericVector & mean, const Eigen::MatrixXd & cholesky);

// evaluate probability density function of a multivariate normal distribution
NumericVector dmvnorm(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance);
NumericVector dmvnorm_cholesky_inverse(const NumericMatrix & x, const NumericVector & mean, const Eigen::MatrixXd & cholesky_inverse);

#endif

