#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
using namespace Rcpp;

DemandLogit::DemandLogit(Eigen::MatrixXd mu_, double u_opt_out_) {
  mu = mu_;
  u_opt_out = u_opt_out_;
  N = mu.cols();
  J = mu.rows();
  N_ = 1/(double) N;

  computed = false;
  rooted = false;

  delta.setZero(J);
  shares.setZero(J,1);
  inc_val = 0.0;
  index_it.setZero(J,1);
  gradient.setZero(J);
  hessian.setZero(J,J);
  obj_val = 0.0;

  //index_i.setZero(J,N);
  //denom_i.setZero(N);
}

//'@export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__new(Eigen::MatrixXd mu_, double u_opt_out_) {
  // create pointer to a object and wrap it as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(new DemandLogit(mu_, u_opt_out_), true);
  return ptr;
}

void DemandLogit::compute(Eigen::VectorXd delta_) {
  double out_index = 0.0;
  double max_index = 0.0;
  double denom_temp = 0.0;

  delta = delta_;

  // zero out the values
  shares.setZero(J,1);
  inc_val = 0.0;
  index_it.setZero(J,1);

  for(int i =0; i < N; i++){
    index_it = mu.col(i).array() + delta.array();

    // safety
    max_index =  index_it.maxCoeff();
    index_it -= max_index;
    out_index = std::exp(u_opt_out - max_index);
    //

    index_it = index_it.exp();

    denom_temp = (index_it.sum() + out_index);

    shares.array() += (index_it / denom_temp);
    inc_val += denom_temp;
  }

  shares.array() *= N_;
  inc_val *= N_;
  computed = true;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__compute(SEXP xp, Eigen::VectorXd delta_) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  ptr->compute(delta_);
  //return the result
  return ptr;
}

Eigen::VectorXd DemandLogit::getShares() {
  return shares;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__getShares(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->getShares();
  //return the result
  return Rcpp::wrap(res);
}

double DemandLogit::getIncValue() {
  return inc_val;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandLogit__getIncValue(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandLogit> ptr(xp);
  //invoke the function
  double res = ptr->getIncValue();
  //return the result
  return wrap(res);
}
