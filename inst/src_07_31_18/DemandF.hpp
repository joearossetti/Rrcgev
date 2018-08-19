#ifndef DEMANDF_HPP
#define DEMANDF_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]

template<const int NUM_PRODUCTS, const int NUM_DRAWS>
class DemandLogitF {
public:
  DemandLogitF() {
    J = NUM_PRODUCTS;
    N = NUM_DRAWS;
    N_ = 1 / (double) N;
    u_opt_out = 0.0;
    computed = false;
    rooted = false;
    inc_val = 0.0;
    //obj_val = 0.0;
  };
  void setS_0(Eigen::VectorXd s_0_) {
    int temp_J = s_0_.size();
    if(temp_J == J) {
      s_0 = s_0_;
    }else{
      Rcpp::stop("wrong dimension");
    }
  };
  Eigen::VectorXd getS_0() {
    Eigen::VectorXd temp = s_0;
    return temp;
  };
  void compute(Eigen::Matrix<double, NUM_PRODUCTS, 1> delta, Eigen::Matrix<double, NUM_PRODUCTS, NUM_DRAWS> mu) {
    double out_index = 0.0;
    double max_index = 0.0;
    double denom_temp = 0.0;

    // zero out the values
    shares.setZero();
    inc_val = 0.0;
    index_it.setZero();

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
  };
  Eigen::VectorXd getShares() {
    if(computed == true){
      return shares;
    }else{
      Rcpp::stop("Did not call compute");
      return shares.setZero();
    }
  };
  double getIncValue() {
    if(computed == true){
      return inc_val;
    }else{
      Rcpp::stop("Did not call compute");
      return 0.0;
    }
  };
protected:
  Eigen::Matrix<double, NUM_PRODUCTS, 1> shares;
  double inc_val;
  Eigen::Matrix<double, NUM_PRODUCTS, 1> s_0;
  //double obj_val;
  //Eigen::Matrix<double, NUM_PRODUCTS, 1> gradient;
  //Eigen::Matrix<double, NUM_PRODUCTS, NUM_PRODUCTS> hessian;
private:
  //Eigen::Matrix<double, NUM_PRODUCTS, NUM_PRODUCTS> mu;
  double u_opt_out;
  double N_;
  int J;
  int N;
  //Eigen::Matrix<double, NUM_PRODUCTS, 1> delta;
  Eigen::Array<double, NUM_PRODUCTS, 1> index_it;
  bool computed;
  bool rooted;
};

#endif /* DEMANDF_HPP */
