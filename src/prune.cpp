#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix prune_measure_(const NumericMatrix & df){
  // df is meant to have a first column with repeat numbers
  // second column made of 1 and 0, indicating (1) whether atom is part of MCMC part or bias correction part
  // third columns with weights
  // and remaining columns (4, 5, ...) with atoms
  // sorted in increasing/lexicographic order
  NumericMatrix newdf(df.nrow(), df.ncol());
  std::fill(newdf.begin(), newdf.end(), 0.);
  // keep track of index of row in newdf in which to write
  int index_row = 0;
  // s will serve as a measure of difference between atoms
  double s = 0;
  newdf(0,_) = df(0,_);
  // loop over rows of df
  for (int i=1; i < df.nrow(); i++){
    // do not prune if different repeat (column 0 contains repeat)
    // if same repeat
    if (df(i,0) == df(i-1,0)){
      // sum absolute value of difference over each component of atom
      s = 0;
      for (int j = 3; j < df.ncol(); j++){
        s += std::abs(df(i,j) - df(i-1,j));
      }
      // if this sum is very small then the atoms are considered equal
      if (s < 1e-20){
        // if so, add up the weights to current weight (column 2 contains weights)
        newdf(index_row, 2)  += df(i, 2);
      } else {
        // new atom, add it to newdf
        index_row++;
        newdf(index_row,_) = df(i,_);
      }
    } else {
      // new atom, add it to newdf
      index_row++;
      newdf(index_row,_) = df(i,_);
    }
  }
  return newdf(Range(0,index_row),_);
}
