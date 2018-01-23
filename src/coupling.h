#ifndef _INCL_CPL_
#define _INCL_CPL_
#include <RcppEigen.h>

using namespace Rcpp;

NumericMatrix gaussian_max_couplingC(const NumericVector & mu1,
                                     const NumericVector & mu2,
                                     const NumericMatrix & Sigma1,
                                     const NumericMatrix & Sigma2);

#endif
