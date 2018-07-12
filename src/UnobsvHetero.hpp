#ifndef UNOBSVHETERO_HPP
#define UNOBSVHETERO_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
//#include <stan/math.hpp>


// [[Rcpp::depends(RcppEigen)]]

class UnobsvHetero {
public:
  UnobsvHetero(Rcpp::NumericMatrix rand_mat,
               Rcpp::LogicalVector which_nln_parmtrs,
               int num_draws, int dim_sigma);
  Rcpp::NumericMatrix getNln_parmtrs_mat(Rcpp::NumericVector short_pars_vec);
  Eigen::MatrixXd getNui_mat(Eigen::MatrixXd sigma_mat);
private:
  int dim_sigma;
  int num_draws;
  Rcpp::LogicalVector which_nln_parmtrs;
  Rcpp::NumericMatrix rand_mat;
};

#endif /* UNOBSVHETERO_HPP */



