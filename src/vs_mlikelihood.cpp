#include <RcppEigen.h>
#include "mvnorm.h"
using namespace Rcpp;

// // [[Rcpp::export]]
// double marginal_likelihood_c_(Eigen::VectorXd selection, const Eigen::MatrixXd & X, const Eigen::VectorXd & Y, double Y2, double g){
//   double l = 0.;
//   int n = X.rows();
//   int p = X.cols();
//   int s = selection.sum();
//   if (s > 0){
//     Eigen::MatrixXd Xselected(n,s);
//     int counter = 0;
//     for (int column = 0; column < p; column++){
//       if (selection(column)){
//         Xselected.col(counter) = X.col(column);
//         counter++;
//       }
//     }
//     l = Y.transpose() * Xselected * ((Xselected.transpose() * Xselected).inverse()) * Xselected.transpose() * Y;
//   } else {
//     l = 0.;
//   }
//   l = -((double) s + 1.) / 2. * log(g + 1.) - (double) n / 2. * log(Y2 - g / (g + 1.) * l);
//   return l;
// }

// [[Rcpp::export]]
double marginal_likelihood_c_2(Eigen::VectorXf selection, const Eigen::MatrixXf & X, const Eigen::VectorXf & Y, double Y2, double g){
  double l = 0.;
  int n = X.rows();
  int p = X.cols();
  int s = selection.sum();
  Eigen::MatrixXf temp;
  if (s > 0){
    Eigen::MatrixXf Xselected(n,s);
    int counter = 0;
    for (int column = 0; column < p; column++){
      if (selection(column)){
        Xselected.col(counter) = X.col(column);
        counter++;
      }
    }
    temp = Y.transpose() * Xselected;
    temp = temp * ((Xselected.transpose() * Xselected).inverse()) * temp.transpose();
    l = temp(0,0);
  } else {
    l = 0.;
  }
  l = -((double) s + 1.) / 2. * log(g + 1.) - (double) n / 2. * log(Y2 - g / (g + 1.) * l);
  return l;
}

// // [[Rcpp::export]]
// double marginal_likelihood_c_3(Eigen::VectorXf selection, const Eigen::MatrixXf & X, const Eigen::VectorXf & Y, double Y2, double g){
//   double l = 0.;
//   int n = X.rows();
//   int p = X.cols();
//   int s = selection.sum();
//   if (s > 0){
//     Eigen::MatrixXf Xselected(n,s);
//     int counter = 0;
//     for (int column = 0; column < p; column++){
//       if (selection(column)){
//         Xselected.col(counter) = X.col(column);
//         counter++;
//       }
//     }
//     Eigen::BDCSVD<Eigen::MatrixXf > svd(Xselected, Eigen::ComputeThinU);
//     Eigen::VectorXf z = svd.matrixU().transpose() * Y.head(s);
//     l = z.dot(z);
//   } else {
//     l = 0.;
//   }
//   l = -((double) s + 1.) / 2. * log(g + 1.) - (double) n / 2. * log(Y2 - g / (g + 1.) * l);
//   return l;
// }


