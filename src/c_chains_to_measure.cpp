#include <Rcpp.h>
using namespace Rcpp;

// This function takes coupled chains, containing chains X (samples1) and Y (samples2)
// and returns a signed empirical measure with atoms and weights.
// lag is an integer >= 1.

// The chains are generated as follows:
// X_0 ~ pi_0, then X_t | X_{t-1} ~ P(X_{t-1}, dot) for t=1,...,lag,
// Y_0 ~ pi_0.

// Then they are coupled, to produce
// (X_{t},Y_{t-lag})|(X_{t-1},Y_{t-lag-1}) ~ bar{P}((X_{t-1},Y_{t-lag-1}), dot)
// for t = lag+1,...,max(m,tau)
// The meeting time is denoted by tau; it is the smallest t such that X_{t} = Y_{t-lag}.

// Thus upon completion, we have
// X_0,..., X_{max(m,tau)}  (of length 1+max(m, tau))
// Y_0,..., Y_{max(m,tau)-lag} (of length 1+max(m, tau) - lag)

// The signed measure approximation, for a lag of one is given by
// (MCMC):             (m-k+1)^{-1} \sum_{t=k}^m delta_{X_t}
// (bias correction):  (m-k+1)^{-1} \sum_{t=k+1}^{tau-1} min(m-k+1, t-k) {delta_{X_t} - delta_{Y_{t-1}}}

// For a general lag, it is given by
// (MCMC):             (m-k+1)^{-1} \sum_{t=k}^m delta_{X_t}
// (bias correction):  (m-k+1)^{-1} \sum_{t=k+lag}^{tau-1} min(m-k+1, ceil((t-k)/lag)) {delta_{X_t} - delta_{Y_{t-lag}}}

// sanity check: it should be the same if lag = 1

// Thus total number of atoms is
// m-k+1                    from the MCMC approximation
// + 2 * max(0, (tau-1) - (k+lag) + 1) from the bias correction term
//  the max above simplifies to max(0, tau - (k + lag))

// Of course some of the atoms are identical but we can deal with that separately,
// by pruning the generated data set later.

// [[Rcpp::export]]
List c_chains_to_measure_as_list_(const List & c_chains, int k, int m){
  int meetingtime = c_chains["meetingtime"];
  NumericMatrix samples1 = c_chains["samples1"];
  NumericMatrix samples2 = c_chains["samples2"];
  int lag = samples1.nrow() - samples2.nrow();
  int size;
  if (((meetingtime-1) - (k+lag) + 1) > 0){
    size = (m - k + 1) + 2 * ((meetingtime-1) - (k+lag) + 1);
  } else {
    size = (m - k + 1);
  }
  // equal weight for each atom in the MCMC part
  double eqweight = 1. / (double) (m - k + 1);
  // number of columns in samples1
  // (should be equal to number of columns in samples2)
  int dimension = samples1.cols();
  // matrix of all atoms
  NumericMatrix atoms(size, dimension);
  // vector of all weights
  NumericVector weights(size);
  // vector indicating which atoms are part of the MCMC measure and which are part of the bias correction
  // TRUE refers to MCMC, FALSE to bias correction
  LogicalVector whichpart(size);
  // fill in MCMC part
  for (int i = 0; i < (m-k+1); i ++){
    atoms(i,_) = samples1(k+i,_);
    weights(i) = eqweight;
    whichpart(i) = true;
  }
  // index keeps track of which row to fill in "atoms"
  int index = m-k+1;
  // now bias correction part
  if (((meetingtime-1) - (k+lag) + 1) > 0){
    for (int time = k+lag; time <= meetingtime-1; time ++){
      atoms(index,_) = samples1(time,_);
      atoms(index+1,_) = samples2(time-lag,_);
      weights(index) = ceil(((double) time - k)/((double) lag));
      if ((m - k + 1) < weights(index)){
        weights(index) = m - k + 1;
      }
      weights(index) = weights(index) * eqweight;
      weights(index+1) = - weights(index);
      whichpart(index) = false;
      whichpart(index+1) = false;
      index += 2;
    }
  }
  return List::create(Named("atoms") = atoms, Named("weights") = weights, Named("MCMC") = whichpart);
}

