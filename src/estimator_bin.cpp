#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// given c_chains, a list produced by the function 'coupled_chains',
// returns estimator of probability of component being between lower and upper
// [[Rcpp::export]]
double estimator_bin(List c_chains, int component, double lower, double upper, int k, int K){
  int meetingtime = c_chains["meetingtime"];
  int iteration = c_chains["iteration"];
  NumericMatrix samples1 = c_chains["samples1"];
  NumericMatrix samples2 = c_chains["samples2"];
  double estimator = 0;
  for (int isample = k; isample <= K; isample ++){
    if (samples1(isample,component-1) > lower && samples1(isample,component-1) < upper){
      estimator += 1;
    }
  }
  if (meetingtime > k + 1){
    double coefficient = 0.;
    double increment = 0;
    for (int isample = k; isample <= std::min(iteration-1, meetingtime-1); isample ++){
      increment = 0;
      coefficient = std::min(isample - k + 1, K - k + 1);
      if (samples1(isample+1,component-1) > lower && samples1(isample+1,component-1) < upper){
        increment += coefficient;
      }
      if (samples2(isample,component-1) > lower && samples2(isample,component-1) < upper){
        increment -= coefficient;
      }
      estimator += increment;
    }
  }
  return estimator / (K - k + 1);
}
