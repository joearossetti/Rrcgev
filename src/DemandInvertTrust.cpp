#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
#include "DemandRoot.hpp"
#include "DemandInvertTrust.hpp"
//#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]


// will need one exported constructor for each type!
//'@export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertTrust_Logit__new(Eigen::VectorXd delta_init, Eigen::VectorXd s_0, Eigen::MatrixXd mu_, double u_opt_out_,
                                             double max_radius, Eigen::VectorXd thresh, Eigen::VectorXd scale,
                                             double tol, int max_count) {
  // create pointer to a object and wrap it as an external pointer
  Rcpp::XPtr< DemandInvertTrust > ptr(new DemandInvertTrust(delta_init, s_0, mu_, u_opt_out_,
                                                       max_radius, thresh, scale, tol, max_count), true);
  return ptr;
}

Eigen::VectorXd DemandInvertTrust::root_cauchy() {
  error = 1000.0;
  int count = 0;
  delta_p = delta;
  int J = delta.size();
  p.setZero(J);
  double rho = 1000.0;
  double value_old = 0.0;

  Rcpp::Rcout << "I started" << std::endl;

  double l1_norm_s_0 = s_0.array().abs().maxCoeff();

  while(true){
    // compute_obj for current delta
    this->compute(delta);
    this->calcValue();
    //this->calcGradient();
    this->calcHessian();
    //save the value of the obj under the old/current delta
    value_old = this->getValue();

    Rcpp::Rcout << "Looping" << std::endl;
    // solve trust region sub-problem
    solve_cauchy();

    // // compute_obj for proposed delta to get rho and error
     this->compute(delta_p);
     this->calcValue();
    // this->calcGradient();
    //
    // rho = (value_old - value) / (value_old - sub_obj_val);
    // error = gradient.array().abs().maxCoeff() / l1_norm_s_0;
    //
    // if(error<tol){
    //   break;
    // }else{
    //   // radius logic
    //   if(rho <= thresh(2)){
    //     // bad approx shrink radius
    //     radius = scale(1) * radius;
    //   }else{
    //     if(rho >= thresh(3) && p.norm()==radius){
    //       // good approx and full step: inc. radius
    //       radius = std::min(scale(2) * radius, max_radius);
    //     }else{
    //       // radius does not change
    //     }
    //   }
    //   // delta iteration logic
    //   if(rho >= thresh(1)){
    //     // accept step
    //     delta = delta_p;
    //   }else{
    //     // bad step resolve with smaller radius (because thresh(1) <= thresh(2))
    //     // delta does not change
    //   }
    // }
    count++;
    if(count>=max_count){
      break;
    }
  }
  Rcpp::Rcout << "The count is " << count << std::endl;
  return delta;
}

//' @export
// [[Rcpp::export]]
RcppExport SEXP DemandInvertTrust_Logit__root_cauchy(SEXP xp) {
  // grap the object as an external pointer
  Rcpp::XPtr< DemandInvertTrust > ptr(xp);
  //invoke the function
  Eigen::VectorXd res = ptr->root_cauchy();
  //return the result
  return Rcpp::wrap(res);
}

void DemandInvertTrust::solve_cauchy() {
  //double gBg = gradient.transpose().dot(hessian * gradient);
  //double norm_g = gradient.norm();
  //Eigen::VectorXd pkc;

  // if(gBg <= 0.0){
  //   p = -(radius / norm_g) * gradient;
  // }else{
  //   if((std::pow(norm_g, 3) / (radius * gBg)) < 1){
  //     p = -(std::pow(norm_g, 2) / (gBg)) * gradient;
  //   }else{
  //     p = -(radius / norm_g) * gradient;
  //   }
  // }
  //delta_p = delta + p;
  //sub_obj_val = p.transpose() * hessian * p + p.dot(gradient) + value;
  delta = delta;
}

