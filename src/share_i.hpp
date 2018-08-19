#ifndef SHARE_I_HPP
#define SHARE_I_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;
using Eigen::DenseBase;
using Eigen::Ref;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
* individual shares are logit
*/
class Logit_i  {
public:
  Logit_i() : u_opt_out(0.0) {};
  Logit_i(double u_opt_out_) : u_opt_out(u_opt_out_) {};
  double u_opt_out;
  void inc_val(const Ref<const VectorXd> u, double& iv) {
    VectorXd y = u.array().exp();
    iv = y.sum() + std::exp(u_opt_out);
    iv = std::log(iv);
  };
  void share(const Ref<const VectorXd> u, double& iv, Ref<VectorXd> s) {
    VectorXd y = u.array().exp();
    iv = y.sum() + std::exp(u_opt_out);
    s = y / iv;
    iv = std::log(iv);
  };
  void jacobian(const Ref<const VectorXd> u, double& iv, Ref<VectorXd> s, Ref<MatrixXd> J) {
    VectorXd y = u.array().exp();
    iv = y.sum() + std::exp(u_opt_out);
    s = y / iv;
    iv = std::log(iv);
    J = s * s.transpose();
    J.diagonal().array() = s.array() * (1 - s.array());
  };
};

/*
 * individual shares are logit
 */
class CNLogit_i  {
public:
  CNLogit_i(MatrixXd A_) : u_opt_out(0.0), A(A_) {};
  CNLogit_i(double u_opt_out_, MatrixXd A_) : u_opt_out(u_opt_out_), A(A_) {};
  double u_opt_out;
  MatrixXd A;
  void inc_val(const Ref<const VectorXd> u, const Ref<const VectorXd> id, const Ref<const VectorXd> mu,
               double& iv) {
    int num_groups = (int) mu.size();
    int num_prods = (int) id.size();
    VectorXd alpha_y_mu;
    VectorXd H;
    H.setZero(num_prods);
    VectorXd y = u.array().exp();
    double weight_mu;

    for(int m = 0; m < num_groups; m++) {
      alpha_y_mu = A.row(m).array() * y.array().pow(mu(m));
      weight_mu = alpha_y_mu.sum();
      weight_mu = std::pow(weight_mu, mu(m)) / weight_mu;
      H.array() += alpha_y_mu.array() * weight_mu;
    }

    iv = H.sum() + std::exp(u_opt_out);
    iv = std::log(iv);
  };
  //void share(const Ref<const VectorXd> u, const Ref<const VectorXd> id, const Ref<const VectorXd> mu,
  //           double& iv, Ref<VectorXd> s) {
  //  VectorXd y = u.array().exp();
  //  iv = y.sum() + std::exp(u_opt_out);
  //  s = y / iv;
  //  iv = std::log(iv);
  //};
  //void jacobian(const Ref<const VectorXd> u, const Ref<const VectorXd> id, const Ref<const VectorXd> mu,
  //              double& iv, Ref<VectorXd> s, Ref<MatrixXd> J) {
  //  VectorXd y = u.array().exp();
  //  iv = y.sum() + std::exp(u_opt_out);
  //  s = y / iv;
  //  iv = std::log(iv);
  // J = s * s.transpose();
  //  J.diagonal().array() = s.array() * (1 - s.array());
 // };
};


#endif /* SHARE_I_HPP */




