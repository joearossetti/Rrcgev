#include <Rcpp.h>
#include <RcppEigen.h>
#include "DemandInvertLogit.hpp"
using namespace Rcpp;

Eigen::VectorXd DemandInvertLogit::invert_fp(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  double error = 1000.0;
  Eigen::ArrayXd delta = delta_0.array();
  Eigen::ArrayXd log_s_diff;

  log_s_diff.setConstant(J, error);
  int count = 0;
  while(true){
    // do update
    log_s_diff = shares_0.array().log() - shares(delta).array().log();
    error = log_s_diff.abs().maxCoeff() / delta.abs().maxCoeff();
    if(error<tol){
      break;
    }else{
      delta += log_s_diff;
    }
    count++;
    if(count>=max_count){
      break;
    }
  }
  Rcout << "The count is " << count << std::endl;
  return delta.matrix();
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertLogit__invert_fp(SEXP xp, Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandInvertLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->invert_fp(shares_0, delta_0, tol, max_count);
  //return the result
  return Rcpp::wrap(res);
}

Eigen::VectorXd DemandInvertLogit::invert_nr(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  double error = 1000.0;
  Eigen::VectorXd delta = delta_0;
  Eigen::VectorXd share_diff;
  Eigen::VectorXd shares_delta;
  Eigen::VectorXd delta_diff;
  Eigen::MatrixXd share_jacobian;

  shares_delta.setZero(J);
  delta_diff.setConstant(J,error);
  share_jacobian.setZero(J,J);
  share_diff.setConstant(J,error);

  double share_zero_l1norm = shares_0.array().abs().maxCoeff();

  int count = 0;
  while(true){

    shares_delta = shares(delta);
    // compute Jacobian
    share_jacobian = shares_delta * shares_delta.transpose();
    share_jacobian.diagonal().array() *= -1.0;
    share_jacobian += shares_delta.asDiagonal();

    // do update
    share_diff = shares_0-shares_delta;
    // solving for roots of shares_delta - shares_0 requires negative function
    delta_diff = share_jacobian.ldlt().solve(share_diff);
    //error = share_diff.array().abs().maxCoeff() / share_zero_l1norm;
    error = delta_diff.array().abs().maxCoeff() / delta.array().abs().maxCoeff();

    if(error<tol){
      break;
    }else{
      delta += delta_diff;
    }
    count++;
    if(count>=max_count){
      break;
    }
  }
  Rcout << "The count is " << count << std::endl;
  return delta;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertLogit__invert_nr(SEXP xp, Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandInvertLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->invert_nr(shares_0, delta_0, tol, max_count);
  //return the result
  return Rcpp::wrap(res);
}

Eigen::VectorXd DemandInvertLogit::invert_nrb(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  double error = 1000.0;
  Eigen::VectorXd delta = delta_0;
  Eigen::VectorXd delta_new;
  Eigen::VectorXd share_diff;
  Eigen::VectorXd share_diff_new;
  Eigen::VectorXd shares_delta;
  Eigen::VectorXd shares_delta_new;
  Eigen::VectorXd delta_diff;
  Eigen::MatrixXd share_jacobian;

  shares_delta.setZero(J);
  shares_delta_new.setZero(J);
  delta_new.setZero(J);
  delta_diff.setConstant(J,error);
  share_jacobian.setZero(J,J);
  share_diff.setConstant(J,error);
  share_diff_new.setConstant(J,error);

  double share_zero_l1norm = shares_0.array().abs().maxCoeff();

  int count = 0;
  while(true){
    // compute chares and see how close to a solution we are
    shares_delta = shares(delta);
    share_diff = (shares_0-shares_delta).array().abs();
    error = share_diff.maxCoeff() / share_zero_l1norm;
    if(error<tol){
      // terminate if near solution
      break;
    }else{
      // do a NR iteration

      // compute Jacobian
      share_jacobian = shares_delta * shares_delta.transpose();
      share_jacobian.diagonal().array() *= -1.0;
      share_jacobian += shares_delta.asDiagonal();

      // do update
      delta_diff = share_jacobian.ldlt().solve(share_diff);

      // check quality of update and adjust
      delta_new = delta + delta_diff;
      shares_delta_new = shares(delta_new);
      share_diff_new = (shares_0-shares_delta_new).array().abs();


      for(int i = 0; i < J; i++){
        if(share_diff_new.coeff(i) > share_diff(i)){
          delta_diff(i) =  std::log(shares_0(i)) - std::log(shares_delta(i));
        }
      }


      delta += delta_diff;
    }
    count++;
    if(count>=max_count){
      break;
    }
  }
  Rcout << "The count is " << count << std::endl;
  return delta;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertLogit__invert_nrb(SEXP xp, Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandInvertLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->invert_nrb(shares_0, delta_0, tol, max_count);
  //return the result
  return Rcpp::wrap(res);
}

Eigen::VectorXd DemandInvertLogit::invert_nrl(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  double error = 1000.0;
  Eigen::VectorXd delta = delta_0;
  Eigen::VectorXd share_diff;
  Eigen::VectorXd shares_delta;
  Eigen::VectorXd delta_diff;
  Eigen::MatrixXd share_jacobian;

  shares_delta.setZero(J);
  delta_diff.setConstant(J,error);
  share_jacobian.setZero(J,J);
  share_diff.setConstant(J,error);

  //double share_zero_l1norm = shares_0.array().abs().maxCoeff();
  double share_zero_l1norm = shares_0.array().log().abs().maxCoeff();

  int count = 0;
  while(true){

    shares_delta = shares(delta);
    // compute Jacobian
    share_jacobian = shares_delta * shares_delta.transpose();
    share_jacobian.diagonal().array() *= -1.0;
    share_jacobian += shares_delta.asDiagonal();

    // switch to log derivative of s0 - shares_delta
    share_jacobian.array().colwise() /= shares_delta.array();

    // do update
    share_diff = shares_0.array().log() - shares_delta.array().log();
    delta_diff = share_jacobian.householderQr().solve(share_diff);
    error = share_diff.array().abs().maxCoeff() / share_zero_l1norm;

    if(error<tol){
      break;
    }else{
      delta += delta_diff;
    }
    count++;
    if(count>=max_count){
      break;
    }
  }
  Rcout << "The count is " << count << std::endl;

  return delta;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertLogit__invert_nrl(SEXP xp, Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandInvertLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->invert_nrl(shares_0, delta_0, tol, max_count);
  //return the result
  return Rcpp::wrap(res);
}

Eigen::VectorXd DemandInvertLogit::invert_fpj(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  double error = 1000.0;
  Eigen::ArrayXd delta = delta_0.array();
  Eigen::ArrayXd log_s_diff;
  Eigen::ArrayXd shares_delta;

  log_s_diff.setConstant(J, error);
  shares_delta.setZero(J);
  int count = 0;
  while(true){
    // do update
    shares_delta = shares(delta).array();
    log_s_diff = (1-shares_delta)*(shares_0.array().log() - shares_delta.log());
    error = log_s_diff.abs().maxCoeff() / delta.abs().maxCoeff();
    if(error<tol){
      break;
    }else{
      delta += log_s_diff;
    }
    count++;
    if(count>=max_count){
      break;
    }
  }
  Rcout << "The count is " << count << std::endl;
  return delta.matrix();
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertLogit__invert_fpj(SEXP xp, Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandInvertLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->invert_fpj(shares_0, delta_0, tol, max_count);
  //return the result
  return Rcpp::wrap(res);
}

Eigen::VectorXd DemandInvertLogit::invert_trcp(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  double error = 1000.0;
  double radius = 1.0;
  double max_radius = 2.0;
  double thresh1 = 0.125;
  double thresh2 = 0.25;
  double thresh3 = 0.75;
  double scale1 = 0.25;
  double scale2 = 2.0;
  Eigen::ArrayXd delta = delta_0.array();
  Eigen::ArrayXd log_s_diff;
  Eigen::ArrayXd shares_delta;

  log_s_diff.setConstant(J, error);
  shares_delta.setZero(J);
  int count = 0;
  while(true){
    // do update
    shares_delta = shares(delta).array();
    log_s_diff = (1-shares_delta)*(shares_0.array().log() - shares_delta.log());
    error = log_s_diff.abs().maxCoeff() / delta.abs().maxCoeff();

    if(error<tol){
      break;
    }else{
      delta += log_s_diff;

    }
    count++;
    if(count>=max_count){
      break;
    }
  }
  Rcout << "The count is " << count << std::endl;
  return delta.matrix();
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertLogit__invert_trcp(SEXP xp, Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count) {
  // grap the object as an external pointer
  Rcpp::XPtr<DemandInvertLogit> ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->invert_trcp(shares_0, delta_0, tol, max_count);
  //return the result
  return Rcpp::wrap(res);
}

