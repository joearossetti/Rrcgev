#ifndef DEMANDROOT_HPP
#define DEMANDROOT_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]


class DemandLogitRoot : public DemandLogit {
public:
  DemandLogitRoot(Eigen::VectorXd s_0, Eigen::MatrixXd mu_, double u_opt_out_) :
  DemandLogit(mu_, u_opt_out_),
  s_0(s_0) {
      log_s_0 = s_0.array().log();
    }
  //
  void compute_obj(Eigen::VectorXd delta);
  Eigen::VectorXd log_diff(Eigen::VectorXd delta);
  double getValue();
  Eigen::VectorXd getGradient();
  Eigen::MatrixXd getHessian();
  void calcValue();
  void calcGradient();
  void calcHessian();
protected:
  double value;
  Eigen::VectorXd gradient;
  Eigen::MatrixXd hessian;
  Eigen::VectorXd s_0;
  Eigen::VectorXd log_s_0;
private:
  Eigen::VectorXd s;
  Eigen::VectorXd delta_;
  double inc_value;
};

#endif /* DEMANDROOT_HPP */
