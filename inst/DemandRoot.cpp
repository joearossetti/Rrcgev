#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
#include "DemandRoot.hpp"
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]

//'@export
// [[Rcpp::export]]
RcppExport SEXP DemandLogitRoot__new(Eigen::VectorXd s_0, Eigen::MatrixXd mu_, double u_opt_out_) {
  // create pointer to a object and wrap it as an external pointer
  Rcpp::XPtr<DemandLogitRoot> ptr(new DemandLogitRoot(s_0, mu_, u_opt_out_), true);
  return ptr;
}

void DemandLogitRoot::compute_obj(Eigen::VectorXd delta) {
  this->compute(delta);
  inc_value = this->getIncValue();
  s = this->getShares();
  delta_ = delta;
  //gradient = s - s_0;
  //hessian = s * s.transpose();
  //hessian.diagonal().array() = s.array() * (1 - s.array());
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogitRoot__compute_obj(SEXP xp, Eigen::VectorXd delta_) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogitRoot> ptr(xp);
  //invoke the function
  ptr->compute_obj(delta_);
  //return the result
  return ptr;
}

void DemandLogitRoot::calcValue() {
  value = inc_value - s_0.dot(delta_);
}

double DemandLogitRoot::getValue() {
  return value;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogitRoot__getValue(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogitRoot> ptr(xp);
  //invoke the function
  ptr->calcValue();
  double res = ptr->getValue();
  //return the result
  return Rcpp::wrap(res);
}

void DemandLogitRoot::calcGradient() {
  gradient = s - s_0;
}

Eigen::VectorXd DemandLogitRoot::getGradient() {
  return gradient;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogitRoot__getGradient(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogitRoot> ptr(xp);
  //invoke the function
  ptr->calcGradient();
  Eigen::VectorXd res = ptr->getGradient();
  //return the result
  return Rcpp::wrap(res);
}

void DemandLogitRoot::calcHessian() {
  hessian = s * s.transpose();
  hessian.diagonal().array() = s.array() * (1 - s.array());
}

Eigen::MatrixXd DemandLogitRoot::getHessian() {
  return hessian;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogitRoot__getHessian(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogitRoot> ptr(xp);
  //invoke the function
  ptr->calcHessian();
  Eigen::MatrixXd res = ptr->getHessian();
  //return the result
  return Rcpp::wrap(res);
}

Eigen::VectorXd DemandLogitRoot::log_diff(Eigen::VectorXd delta) {
  this->compute(delta);
  value = this->getIncValue();
  s = this->getShares();
  delta_ = delta;
  return (s.array().log() - log_s_0.array()).matrix();
}

