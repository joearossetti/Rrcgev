#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]


void times_two(Eigen::Ref<Eigen::VectorXd> a) {
  a = a*2;
};

//' @export
// [[Rcpp::export]]
Eigen::VectorXd vec_sub(Eigen::VectorXd x, NumericVector start_vec, NumericVector end_vec, int M) {
  int temp_length = 0;
  for(int m=0; m < M; m++) {
    temp_length = end_vec[m] - start_vec[m] + 1;
    //x.segment(start_vec[m], temp_length) *= 2;
    times_two( x.segment(start_vec[m], temp_length) );
  }

  return x;
}

