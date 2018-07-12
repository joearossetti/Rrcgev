#ifndef DEMAND_HPP
#define DEMAND_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]

class DemandLogit {
public:
  DemandLogit(Eigen::MatrixXd mu_, double u_opt_out_);
  void compute(Eigen::VectorXd delta_);
  Eigen::VectorXd getShares();
  double getIncValue();
  void setS_0(Eigen::VectorXd s_0);
  void calcObj_val();
  void calcGradient();
  void calcHessian();
  Eigen::VectorXd getGradient();
  Eigen::MatrixXd getHessian();
  double getObj_val();
protected:
  Eigen::VectorXd shares;
  double inc_val;
  Eigen::VectorXd s_0;
  double obj_val;
  Eigen::VectorXd gradient;
  Eigen::MatrixXd hessian;
private:
  Eigen::VectorXd delta;
  double u_opt_out;
  double N_;
  int J;
  int N;
  Eigen::MatrixXd mu;
  Eigen::ArrayXd index_it;
  bool computed;
  bool rooted;
};

#endif /* DEMAND_HPP */



