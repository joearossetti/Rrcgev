#include <Rcpp.h>
//#include <stan/math.hpp>
#include "UnobsvHetero.hpp"
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

UnobsvHetero::UnobsvHetero(Rcpp::NumericMatrix rand_mat_,
                           Rcpp::LogicalVector which_nln_parmtrs_,
                           int num_draws_, int dim_sigma_) {

  num_draws = num_draws_;
  dim_sigma = dim_sigma_;
  which_nln_parmtrs = which_nln_parmtrs_;
  rand_mat = rand_mat_;
}

// external pointer to uniform object
//' @export
// [[Rcpp::export]]
RcppExport SEXP UnobsvHetero__new(Rcpp::NumericMatrix rand_mat_,
                                  Rcpp::LogicalVector which_nln_parmtrs_,
                                  int num_draws_, int dim_sigma_) {

  // create pointer to a uniform object and wrap it as an external pointer
  Rcpp::XPtr<UnobsvHetero> ptr(new UnobsvHetero(rand_mat_, which_nln_parmtrs_,num_draws_,dim_sigma_), true);

  return ptr;
}

// method to return the nln_mat
Rcpp::NumericMatrix UnobsvHetero::getNln_parmtrs_mat(Rcpp::NumericVector short_pars_vec) {
  Rcpp::NumericMatrix nln_mat(dim_sigma, dim_sigma);

  int index = 0;
  int j = 0;
  for(NumericMatrix::iterator iter = nln_mat.begin(); iter != nln_mat.end(); iter++){
    if(which_nln_parmtrs[index]==true){
      *iter = short_pars_vec[j];
      j++;
    }else{
      *iter = 0.0;
    }
    index++;
  }

  return nln_mat;
}

// invoke the nln_mat method
//' @export
// [[Rcpp::export]]
RcppExport SEXP UnobsvHetero__getNln_parmtrs_mat(SEXP xp, Rcpp::NumericVector short_pars_vec_) {
  // grap the object as an external pointer
  Rcpp::XPtr<UnobsvHetero> ptr(xp);


  //invoke the function
  NumericMatrix res = ptr->getNln_parmtrs_mat(short_pars_vec_);

  //return the result
  return res;
}


// method to get the nui_mat method
Eigen::MatrixXd UnobsvHetero::getNui_mat(Eigen::MatrixXd sigma_mat) {
  Eigen::LLT<Eigen::MatrixXd> lltOf_sigma_mat(sigma_mat);

  return lltOf_sigma_mat.matrixL() * (Rcpp::as<Eigen::MatrixXd>(rand_mat));
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP UnobsvHetero__getNui_mat(SEXP xp, Eigen::MatrixXd sigma_mat_) {
  // grap the object as an external pointer
  Rcpp::XPtr<UnobsvHetero> ptr(xp);
  //invoke the function
  Eigen::MatrixXd res = ptr->getNui_mat(sigma_mat_);
  //return the result
  return Rcpp::wrap(res);
}
