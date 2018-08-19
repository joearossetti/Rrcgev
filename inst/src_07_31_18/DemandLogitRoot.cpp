#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

void DemandLogit::setS_0(Eigen::VectorXd s_0_) {
  s_0 = s_0_;
  rooted = true;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__setS_0(SEXP xp, Eigen::VectorXd s_0_) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  ptr->setS_0(s_0_);
  //return the result
  return ptr;
}

void DemandLogit::calcObj_val() {
  if(rooted == true && computed == true){
    obj_val = inc_val - delta.dot(s_0);
  }else{
    Rcpp::Rcout << "not computed or not rooted" << std::endl;
  }
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__calcObj_val(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  ptr->calcObj_val();
  //return the result
  return ptr;
}

void DemandLogit::calcGradient() {
  if(rooted == true && computed == true){
    gradient = shares - s_0;
  }else{
    Rcpp::Rcout << "not computed or not rooted" << std::endl;
  }
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__calcGradient(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  ptr->calcGradient();
  //return the result
  return ptr;
}

void DemandLogit::calcHessian() {
  if(rooted == true && computed == true){
    hessian = shares * shares.transpose();
    hessian.diagonal().array() = shares.array() * (1 - shares.array());
  }else{
    Rcpp::Rcout << "not computed or not rooted" << std::endl;
  }
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__calcHessian(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  ptr->calcHessian();
  //return the result
  return ptr;
}

Eigen::VectorXd DemandLogit::getGradient() {
  if(rooted == true && computed == true){
    return gradient;
  }else{
    Rcpp::Rcout << "not computed or not rooted" << std::endl;
    return gradient;
  }
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__getGradient(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->getGradient();
  //return the result
  return wrap(res);
}

Eigen::MatrixXd DemandLogit::getHessian() {
  if(rooted == true && computed == true){
    return hessian;
  }else{
    Rcpp::Rcout << "not computed or not rooted" << std::endl;
    return hessian;
  }
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__getHessian(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  Eigen::MatrixXd res = ptr->getHessian();
  //return the result
  return wrap(res);
}

double DemandLogit::getObj_val() {
  if(rooted == true && computed == true){
    return obj_val;
  }else{
    Rcpp::Rcout << "not computed or not rooted" << std::endl;
    return obj_val;
  }
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__getObj_val(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  double res = ptr->getObj_val();
  //return the result
  return wrap(res);
}
