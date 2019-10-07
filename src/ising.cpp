#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;

// given a square grid x = {x_i}, e.g. 32x32 or something
// with each x_i in {-1,+1}
// the following computes sum_{i~j} x_i x_j
// where i~j denotes the neighboring relation; here we use perodic boundary condition
// so (0,0) is neighbor with (32,0) in a 32x32 grid, for instance
// [[Rcpp::export]]
int ising_sum_(const IntegerMatrix & state){
  int size = state.rows();
  int s = 0;
  // below, (i1,j1) will index a neighbor of (i,j)
  int i1;
  int j1;
  for (int i = 0; i < size; i++){
    if (i == (size - 1)){
      i1 = 0;
    } else {
      i1 = i+1;
    }
    for (int j = 0; j < size; j++){
      if (j == (size - 1)){
        j1=0;
      } else {
        j1 = j+1;
      }
      // sums over two neighbors ("top right" neighbors
      // of current location, if counting from bottom left corner)
      s += state(i,j) * (state(i,j1) + state(i1,j));
    }
  }
  return s;
}

// given a square grid 'state'
// and a vector proba_beta of values between 0 and 1
// corresponding to probabilities of drawing a +1 given the neighbors sum to
// -4, -2, 0, 2, 4 (respectively, so proba_beta should contain 5 values)
// then the following functions perform a sweep over all locations in the grid and performs a flip
// [[Rcpp::export]]
IntegerMatrix ising_gibbs_sweep_(IntegerMatrix state, NumericVector proba_beta){
  RNGScope scope;
  int size = state.rows();
  int s;
  int itop, ibottom, jright, jleft;
  for (int i = 0; i < size; i++){
    for (int j = 0; j < size; j++){
      s = 0;
      itop = (i+1) % size;
      ibottom = ((i + size - 1) % size);
      jright = (j+1) % size;
      jleft = (j + size - 1) % size;
      s += state(itop, j) + state(ibottom, j) + state(i, jright) + state(i, jleft);
      GetRNGstate();
      state(i,j) = 2*((runif(1))(0) < proba_beta((s+4)/2)) - 1;
      PutRNGstate();
    }
  }
  return state;
}

// coupled version of single-site Gibbs update,
// where the strategy is to maximally couple each conditional update
// [[Rcpp::export]]
List ising_coupled_gibbs_sweep_(IntegerMatrix state1, IntegerMatrix state2, NumericVector proba_beta){
  RNGScope scope;
  int size = state1.rows();
  int s1;
  int s2;
  int itop,ibottom,jright,jleft;
  IntegerVector x(2);
  for (int i = 0; i < size; i++){
    for (int j = 0; j < size; j++){
      s1 = 0;
      s2 = 0;
      itop = (i+1) % size;
      ibottom = ((i + size - 1) % size);
      jright = (j+1) % size;
      jleft = (j + size - 1) % size;
      s1 += state1(itop, j) + state1(ibottom, j) + state1(i, jright) + state1(i, jleft);
      s2 += state2(itop, j) + state2(ibottom, j) + state2(i, jright) + state2(i, jleft);
      GetRNGstate();
      double u_ = (runif(1))(0);
      PutRNGstate();
      state1(i,j) = 2*(u_ < proba_beta((s1+4)/2)) - 1;
      state2(i,j) = 2*(u_ < proba_beta((s2+4)/2)) - 1;
    }
  }
  return List::create(Named("state1") = state1, Named("state2") = state2);
}


