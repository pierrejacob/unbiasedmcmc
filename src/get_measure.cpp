#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List get_measure_(const List & c_chains, int k, int m){
  int meetingtime = c_chains["meetingtime"];
  int size;
  if (meetingtime > k){
    size = (m - k + 1) + 2 * (meetingtime - k);
  } else {
    size = (m - k + 1);
  }
  double eqweight = 1. / (m - k + 1);
  NumericMatrix samples1 = as<NumericMatrix>(c_chains["samples1"]);
  NumericMatrix samples2 = as<NumericMatrix>(c_chains["samples2"]);
  int dim = samples1.cols();
  NumericMatrix atoms(size, dim);
  NumericVector weights(size);
  for (int i = 0; i < (m-k+1); i ++){
    atoms(i,_) = samples1(k+i,_);
    weights(i) = eqweight;
  }
  int index = m-k+1;
  if (meetingtime < k + 1){
    // nothing
  } else {
    int t = k;
    while (t < meetingtime){
      atoms(index,_) = samples1(t + 1,_);
      atoms(index+1,_) = samples2(t,_);
      weights(index) = (t - k + 1.) * eqweight;
      if (weights(index) > 1){
        weights(index) = 1;
      }
      weights(index+1) = - weights(index);
      t ++;
      index += 2;
    }
  }
  return List::create(Named("atoms") = atoms, Named("weights") = weights);
}
