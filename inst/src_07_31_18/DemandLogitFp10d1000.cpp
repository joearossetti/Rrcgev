#include <Rcpp.h>
#include <RcppEigen.h>
#include "DemandF.hpp"
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

const int num_products_ = 10;
const int num_draws_ = 1000;

typedef DemandLogitF<num_products_, num_draws_> Dlfp10d1000;

//'@export
// [[Rcpp::export]]
RcppExport SEXP Dlfp10d1000__new() {
  // create pointer to a object and wrap it as an external pointer
  Rcpp::XPtr<Dlfp10d1000> ptr(new Dlfp10d1000(), true);
  return ptr;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP Dlfp10d1000__setS_0(SEXP xp, Eigen::VectorXd s_0_) {
  // grap the object as an external pointer
  Rcpp::XPtr<Dlfp10d1000> ptr(xp);
  //invoke the function
  ptr->setS_0(s_0_);
  //return the result
  return ptr;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP Dlfp10d1000__getS_0(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<Dlfp10d1000> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->getS_0();
  //return the result
  return Rcpp::wrap(res);
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP Dlfp10d1000__compute(SEXP xp, Eigen::VectorXd delta_, Eigen::MatrixXd mu_) {
  // grap the object as an external pointer
  Rcpp::XPtr<Dlfp10d1000> ptr(xp);

  Eigen::Ref< Eigen::Matrix<double, 10, 1> > delta_f(delta_);
  Eigen::Ref< Eigen::Matrix<double, 10, 1000> > mu_f(mu_);

  //invoke the function
  ptr->compute(delta_f, mu_f);
  //return the result
  return ptr;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP Dlfp10d1000__getShares(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<Dlfp10d1000> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->getShares();
  //return the result
  return Rcpp::wrap(res);
}
