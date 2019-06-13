#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
using namespace std;

// given c_chains, a list produced by the function 'coupled_chains',
// returns estimator of probability of component being between lower and upper
// [[Rcpp::export]]
double estimator_bin_(List c_chains, int component, double lower, double upper, int k, int m, int lag){
  int meetingtime = c_chains["meetingtime"];
  NumericMatrix samples1 = c_chains["samples1"];
  NumericMatrix samples2 = c_chains["samples2"];
  double estimator = 0;
  for (int isample = k; isample <= m; isample ++){
    if (samples1(isample,component-1) > lower && samples1(isample,component-1) < upper){
      estimator += 1;
    }
  }
  // next, add bias correction terms
  if (meetingtime > k + lag){
    double coefficient = 0.;
    double increment = 0.;
    for (int time = k+lag; time <= meetingtime-1; time ++){
      increment = 0.;
      // compute min between m - k + 1 and ceiling ((t - k) / L)
      coefficient = ceil(((double) time - k)/((double) lag));
      if ((m - k + 1) < coefficient){
        coefficient = m - k + 1.;
      }
      if (samples1(time,component-1) > lower && samples1(time,component-1) < upper){
        increment += coefficient;
      }
      if (samples2(time-lag,component-1) > lower && samples2(time-lag,component-1) < upper){
        increment -= coefficient;
      }
      estimator += increment;
    }
  }
  return estimator / (m - k + 1.);
}
