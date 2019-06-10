#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// This function takes a vector of zeros and ones,
// with at least one zero and one one
// and samples uniformly among the indices pointing to zeros,
// and uniformly among the indices pointing to ones, and
// outputs the pair
// [[Rcpp::export]]
IntegerVector sample_pair01(const NumericVector & selection){
  RNGScope scope;
  int l = selection.size();
  int sumones = sum(selection);
  IntegerVector indices(2);
  indices(0) = -1;
  indices(1) = -1;
  GetRNGstate();
  NumericVector us = runif(2);
  PutRNGstate();
  double u0 = us(0);
  double u1 = us(1);
  double w0 = 1. / (l - sumones);
  double w1 = 1. / sumones;
  double csw0 = 0.;
  double csw1 = 0.;
  for (int k = 0; k < l; k++){
    if (selection(k) == 0){
      csw0 += w0;
      if (indices(0) < 0 && csw0 > u0){
        indices(0) = k+1;
      }
    } else {
      csw1 += w1;
      if (indices(1) < 0 && csw1 > u1){
        indices(1) = k+1;
      }
    }
  }
  return indices;
}

