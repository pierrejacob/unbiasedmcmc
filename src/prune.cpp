#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix prune_(const NumericMatrix & df){
  // df is meant to have a first column with repeat numbers
  // second columns with weights
  // and remaining columns (3, 4, ...) with the atoms
  NumericMatrix copieddf = clone(df);
  double s = 0;
  for (int i=1; i < copieddf.rows(); i++){
    // do not prune if different repeat
    if (copieddf(i,0) == copieddf(i-1,0)){
      s = 0;
      // check if successive rows have identical atoms
      for (int j = 2; j < copieddf.cols(); j++){
        s += std::abs(copieddf(i,j) - copieddf(i-1,j));
      }
      if (s < 1e-20){
        // if so, add up the weights
        copieddf(i, 1)  += copieddf(i-1, 1);
        copieddf(i-1, 1) = 0.;
      }
    }
  }
  return copieddf;
}
