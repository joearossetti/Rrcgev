#include <Rcpp.h>
#include <RcppEigen.h>
#include "Demand.hpp"
#include <unsupported/Eigen/NonLinearOptimization>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


struct DemandInvFunctor : public DemandLogit
{
  // constructor
  DemandInvFunctor(Eigen::MatrixXd mu_, double u_opt_out_, Eigen::VectorXd s_0_) :
    DemandLogit(mu_, u_opt_out_) {
      n = mu_.rows();
      m = n;
      DemandLogit::setS_0(s_0_);
  }

  // Compute 'm' errors, one for each data point, for the given parameter values in 'x'
  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
  {
    DemandLogit::compute(x);
    DemandLogit::calcGradient();

    //Rcpp::Rcout << "calc gradient" << std::endl;

    fvec = gradient;

    return 0;
  }

  // Compute the jacobian of the errors
  int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac)
  {
    DemandLogit::compute(x);
    DemandLogit::calcHessian();

    //Rcpp::Rcout << "calc hessian" << std::endl;

    fjac = hessian;

    return 0;
  }

  // Number of data points, i.e. values.
  int m;

  // Returns 'm', the number of values.
  int values() const { return m; }

  // The number of parameters, i.e. inputs.
  int n;

  // Returns 'n', the number of inputs.
  int inputs() const { return n; }

};

//' @export
// [[Rcpp::export]]
NumericVector root_lm(Eigen::VectorXd delta_init, Eigen::MatrixXd mu_, double u_opt_out_, Eigen::VectorXd s_0_) {
  Eigen::VectorXd x = delta_init;
  DemandInvFunctor my_functor(mu_, u_opt_out_, s_0_);

  Eigen::LevenbergMarquardt<DemandInvFunctor, double> lm(my_functor);
  //lm.setEpsilon(my_tol);
  lm.parameters.maxfev = 1000;
  lm.parameters.xtol = 2.0e-16;
  lm.minimize(x);

  //int count = *lm.iterations();

  return wrap(x);
}

