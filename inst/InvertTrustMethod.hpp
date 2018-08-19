#ifndef INVERTTRUSTMETHOD_HPP
#define INVERTTRUSTMETHOD_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]

class InvertTrustMethod : public DemandLogit {
public:
  InvertTrustMethod(Eigen::VectorXd delta_init, double init_radius, Eigen::VectorXd s_0, Eigen::MatrixXd mu_, double u_opt_out_,
                    double max_radius, Eigen::VectorXd thresh, Eigen::VectorXd scale,
                    double tol, int max_count) :
  DemandLogit(mu_, u_opt_out_), max_radius(max_radius), radius(init_radius),
  thresh(thresh), scale(scale), tol(tol), max_count(max_count), delta(delta_init) { setS_0(s_0); };
  Eigen::VectorXd root_cauchy();
  void solve_cauchy();
  //Eigen::VectorXd solve_dogleg();
private:
  double max_radius;
  double tol;
  int max_count;
  Eigen::VectorXd thresh;
  Eigen::VectorXd scale;
  double radius;
  double error;
  double sub_obj_val;
  Eigen::VectorXd delta;
  Eigen::VectorXd delta_p;
  Eigen::VectorXd p;
};


#endif /* INVERTTRUSTMETHOD_HPP */
