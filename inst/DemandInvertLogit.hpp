#ifndef NOBSVHETERO_HPP
#define NOBSVHETERO_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]


class DemandInvertLogit {
public:
  DemandInvertLogit();
  Eigen::VectorXd invert_fp(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count);
  Eigen::VectorXd invert_fpj(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count);
  Eigen::VectorXd invert_nr(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count);
  Eigen::VectorXd invert_nrb(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count);
  Eigen::VectorXd invert_nrl(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count);
  Eigen::VectorXd invert_trcp(Eigen::VectorXd shares_0, Eigen::VectorXd delta_0, double tol, int max_count);
private:
};


#endif /* NOBSVHETERO_HPP */
