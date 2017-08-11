#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix sigma_(const NumericMatrix & X, const NumericVector & w){
  int n = X.rows();
  int p = X.cols();
  NumericMatrix inv_sigma(p,p);
  for (int j1 = 0; j1 < p; j1 ++){
    for (int j2 = j1; j2 < p; j2 ++){
      inv_sigma(j1,j2) = 0;
      for (int i = 0; i < n; i++){
        inv_sigma(j1,j2) = inv_sigma(j1,j2) + X(i,j1) * X(i,j2) * w(i);
      }
    }
  }
  for (int j1 = 1; j1 < p; j1 ++){
    for (int j2 = 0; j2 < j1; j2 ++){
      inv_sigma(j1,j2) = inv_sigma(j2,j1);
    }
  }
  return inv_sigma;
}

// The following function computes m(omega) and Sigma(omega)... (or what we really need instead)
// it returns m (= m(omega)), Sigma_inverse = Sigma(omega)^{-1},
// as well as Cholesky_inverse and Cholesky that are such that
// Cholesky_inverse is the lower triangular matrix L, in the decomposition Sigma^{-1} = L L^T
// whereas Cholesky is the lower triangular matrix Ltilde in the decomposition Sigma = Ltilde^T Ltilde

// [[Rcpp::export]]
List m_sigma_function_(const Eigen::Map<Eigen::MatrixXd>  & omega,
                       const Eigen::Map<Eigen::MatrixXd>  & X,
                       const Eigen::Map<Eigen::MatrixXd>  & invB,
                       const Eigen::Map<Eigen::VectorXd>  & KTkappaplusinvBtimesb){
  int n = X.rows();
  int p = X.cols();
  // The matrix A stores XT Omega X + B^{-1}, that is, Sigma^{-1}
  Eigen::MatrixXd A(p,p);
  for (int j1 = 0; j1 < p; j1 ++){
    for (int j2 = j1; j2 < p; j2 ++){
      A(j1,j2) = invB(j1, j2);
      for (int i = 0; i < n; i++){
        A(j1,j2) = A(j1,j2) + X(i,j1) * X(i,j2) * omega(i);
      }
      A(j2,j1) = A(j1,j2);
    }
  }
  Eigen::LLT<Eigen::MatrixXd> lltofA(A);
  Eigen::MatrixXd lower = lltofA.matrixL();
  Eigen::VectorXd x = lltofA.solve(KTkappaplusinvBtimesb);
  return List::create(Named("m")=x,
                      Named("Sigma_inverse") = A,
                      Named("Cholesky_inverse") = lower,
                      Named("Cholesky") = lower.inverse());
}
