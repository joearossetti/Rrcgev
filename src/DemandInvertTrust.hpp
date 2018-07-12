#ifndef DEMANDINVERTTRUST_HPP
#define DEMANDINVERTTRUST_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
#include "DemandRoot.hpp"
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]

// template <class DemandType>
// class DemandInvertTrust : public DemandType{
// public:
//   DemandInvertTrust();
//   Eigen::VectorXd root_cauchy();
//   void solve_cauchy();
//   //Eigen::VectorXd solve_dogleg();
// private:
//   double max_radius;
//   Eigen::Vector3d thresh;
//   Eigen::Vector2d scale;
//   double tol;
//   int max_count;
//   double radius;
//   double error;
//   Eigen::VectorXd delta;
//   Eigen::VectorXd delta_p;
//   Eigen::VectorXd p;
//   double sub_obj_val;
// };

//template <>
class DemandInvertTrust : public DemandLogitRoot{
public:
  DemandInvertTrust(Eigen::VectorXd delta_init, Eigen::VectorXd s_0, Eigen::MatrixXd mu_, double u_opt_out_,
                    double max_radius, Eigen::VectorXd thresh, Eigen::VectorXd scale,
                    double tol, int max_count) :
  DemandLogitRoot(s_0, mu_, u_opt_out_), max_radius(max_radius),
  thresh(thresh), scale(scale), tol(tol), max_count(max_count), delta(delta_init) {};
  Eigen::VectorXd root_cauchy();
  void solve_cauchy();
  //Eigen::VectorXd solve_dogleg();
private:
  double max_radius;
  Eigen::VectorXd thresh;
  Eigen::VectorXd scale;
  double tol;
  int max_count;
  double radius;
  double error;
  Eigen::VectorXd delta;
  Eigen::VectorXd delta_p;
  Eigen::VectorXd p;
  double sub_obj_val;
};


#endif /* DEMANDINVERTTRUST_HPP */
