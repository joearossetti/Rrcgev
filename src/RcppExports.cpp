// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// logit_i_test_shares
Eigen::VectorXd logit_i_test_shares(Eigen::VectorXd utils);
RcppExport SEXP _Rrcgev_logit_i_test_shares(SEXP utilsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type utils(utilsSEXP);
    rcpp_result_gen = Rcpp::wrap(logit_i_test_shares(utils));
    return rcpp_result_gen;
END_RCPP
}
// logit_i_test_val
double logit_i_test_val(Eigen::VectorXd utils);
RcppExport SEXP _Rrcgev_logit_i_test_val(SEXP utilsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type utils(utilsSEXP);
    rcpp_result_gen = Rcpp::wrap(logit_i_test_val(utils));
    return rcpp_result_gen;
END_RCPP
}
// logit_i_test_jac
Eigen::MatrixXd logit_i_test_jac(Eigen::VectorXd utils);
RcppExport SEXP _Rrcgev_logit_i_test_jac(SEXP utilsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type utils(utilsSEXP);
    rcpp_result_gen = Rcpp::wrap(logit_i_test_jac(utils));
    return rcpp_result_gen;
END_RCPP
}
// vec_sub
Eigen::VectorXd vec_sub(Eigen::VectorXd x, NumericVector start_vec, NumericVector end_vec, int M);
RcppExport SEXP _Rrcgev_vec_sub(SEXP xSEXP, SEXP start_vecSEXP, SEXP end_vecSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type start_vec(start_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type end_vec(end_vecSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(vec_sub(x, start_vec, end_vec, M));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rrcgev_logit_i_test_shares", (DL_FUNC) &_Rrcgev_logit_i_test_shares, 1},
    {"_Rrcgev_logit_i_test_val", (DL_FUNC) &_Rrcgev_logit_i_test_val, 1},
    {"_Rrcgev_logit_i_test_jac", (DL_FUNC) &_Rrcgev_logit_i_test_jac, 1},
    {"_Rrcgev_vec_sub", (DL_FUNC) &_Rrcgev_vec_sub, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rrcgev(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
