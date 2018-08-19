#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "share_i.hpp"
using namespace Rcpp;
using namespace Eigen;
using Eigen::DenseBase;
using Eigen::Ref;
using Eigen::MatrixXd;
using Eigen::VectorXd;

//'@export
// [[Rcpp::export]]
RcppExport SEXP Logit_i__new() {
  // create pointer to a object and wrap it as an external pointer
  Rcpp::XPtr<Logit_i> ptr(new Logit_i(), true);
  return ptr;
}

//'@export
// [[Rcpp::export]]
RcppExport SEXP Logit_i_out__new(double u_opt_out_) {
  // create pointer to a object and wrap it as an external pointer
  Rcpp::XPtr<Logit_i> ptr(new Logit_i(u_opt_out_), true);
  return ptr;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP Logit_i__share(SEXP xp, Eigen::VectorXd u) {
  // grap the object as an external pointer
  Rcpp::XPtr<Logit_i> ptr(xp);

  Eigen::VectorXd s_;

  //invoke the function
  ptr->share(u, s_);
  //return the result
  return Rcpp::wrap(s_);
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP Logit_i__jacobian(SEXP xp, Eigen::VectorXd u) {
  // grap the object as an external pointer
  Rcpp::XPtr<Logit_i> ptr(xp);

  Eigen::VectorXd s_;
  Eigen::MatrixXd J_;

  //invoke the function
  ptr->jacobian(u, s_, J_);
  //return the result
  return Rcpp::wrap(J_);
}
