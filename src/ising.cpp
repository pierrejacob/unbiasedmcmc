#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;

// samples two bernoulli with probabilities p1 and p2
// such that it is a maximal coupling
// and the samples take values in {-1,+1}
// [[Rcpp::export]]
IntegerVector rcouplbern2(double p1, double p2){
  RNGScope scope;
  IntegerVector x(2);
  NumericVector u(1);
  double p10 = 1-p1;
  double p20 = 1-p2;
  double minoriz0, minoriz1;
  double alpha;
  if (p10 < p20){
    minoriz0 = p10;
    minoriz1 = p2;
  } else {
    minoriz0 = p20;
    minoriz1 = p1;
  }
  alpha = minoriz0 + minoriz1;
  u = runif(1);
  if (u(0) < alpha){
    u = runif(1);
    x(0) = (u(0) < minoriz1/alpha);
    x(1) = x(0);
  } else {
    u = runif(1);
    x(0) = (u(0) < (p1-minoriz1)/(1.-alpha));
    u = runif(1);
    x(1) = (u(0) < (p2-minoriz1)/(1.-alpha));
  }
  return((2*x - 1));
}

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
      state(i,j) = 2*((runif(1))(0) < proba_beta((s+4)/2)) - 1;
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
  NumericVector u(1);
  double p1, p10, p2, p20;
  double minoriz0, minoriz1;
  double alpha;
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
      p1 = proba_beta((s1+4)/2);
      p10 = 1-p1;
      p2 = proba_beta((s2+4)/2);
      p20 = 1-p2;
      //
      if (p10 < p20){
        minoriz0 = p10;
        minoriz1 = p2;
      } else {
        minoriz0 = p20;
        minoriz1 = p1;
      }
      alpha = minoriz0 + minoriz1;
      u = runif(1);
      if (u(0) < alpha){
        u = runif(1);
        x(0) = (u(0) < minoriz1/alpha);
        x(1) = x(0);
      } else {
        u = runif(1);
        x(0) = (u(0) < (p1-minoriz1)/(1.-alpha));
        u = runif(1);
        x(1) = (u(0) < (p2-minoriz1)/(1.-alpha));
      }
      state1(i,j) = 2*x(0) - 1;
      state2(i,j) = 2*x(1) - 1;
    }
  }
  return List::create(Named("state1") = state1, Named("state2") = state2);
}
// List coupled_gibbs_sweep(IntegerMatrix state1, IntegerMatrix state2, NumericVector proba_beta){
//   RNGScope scope;
//   int size = state1.rows();
//   int s1;
//   int s2;
//   int itop,ibottom,jright,jleft;
//   for (int i = 0; i < size; i++){
//     for (int j = 0; j < size; j++){
//       s1 = 0;
//       s2 = 0;
//       itop = (i+1) % size;
//       ibottom = ((i + size - 1) % size);
//       jright = (j+1) % size;
//       jleft = (j + size - 1) % size;
//       s1 += state1(itop, j) + state1(ibottom, j) + state1(i, jright) + state1(i, jleft);
//       s2 += state2(itop, j) + state2(ibottom, j) + state2(i, jright) + state2(i, jleft);
//       IntegerVector xx = rcouplbern2(proba_beta((s1+4)/2), proba_beta((s2+4)/2));
//       state1(i,j) = xx(0);
//       state2(i,j) = xx(1);
//     }
//   }
//   return List::create(Named("state1")=state1, Named("state2") = state2);
// }


// [[Rcpp::export]]
int ising_two2one_(int ix, int iy, int size){
  return ix + iy * size;
}

// [[Rcpp::export]]
IntegerVector ising_one2two_(int location, int size){
  IntegerVector is(2);
  is(0) = location % size;
  is(1) = location / size;
  return is;
}

// [[Rcpp::export]]
int ising_locationneighbour_(int location, int ineighbour, int size){
  IntegerVector location2 = ising_one2two_(location, size);
  if (ineighbour == 0){ // neighbour on the right
    location2(0) += 1;
    if (location2(0) == size){
      location2(0) = 0;
    }
  }
  if (ineighbour == 1){ // neighbour on the left
    location2(0) -= 1;
    if (location2(0) == -1){
      location2(0) = size-1;
    }
  }
  if (ineighbour == 2){ // neighbour on the top
    location2(1) += 1;
    if (location2(1) == size){
      location2(1) = 0;
    }
  }
  if (ineighbour == 3){ // neighbour on the bottom
    location2(1) -= 1;
    if (location2(1) == -1){
      location2(1) = size-1;
    }
  }
  return ising_two2one_(location2(0), location2(1), size);
}
