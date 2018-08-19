#include <Rcpp.h>
#include <RcppEigen.h>
#include "share_i.hpp"
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::VectorXd logit_i_test_shares(Eigen::VectorXd utils) {
  Logit_i demand;
  Eigen::VectorXd shares = utils;
  double value = 0.0;

  demand.share(utils, value, shares);

  return shares;
}

// [[Rcpp::export]]
double logit_i_test_val(Eigen::VectorXd utils) {
  Logit_i demand;
  double value = 0.0;

  demand.inc_val(utils, value);

  return value;
}

// [[Rcpp::export]]
Eigen::MatrixXd logit_i_test_jac(Eigen::VectorXd utils) {
  Logit_i demand;
  Eigen::VectorXd shares = utils;
  double value = 0.0;
  Eigen::MatrixXd Jac = shares * shares.transpose();

  demand.jacobian(utils, value, shares, Jac);

  return Jac;
}


