// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// fit_twopart_cpp
Rcpp::List fit_twopart_cpp(const Rcpp::NumericMatrix& X_, const Rcpp::NumericVector& Z_, const Rcpp::NumericMatrix& Xs_, const Rcpp::NumericVector& S_, const Rcpp::IntegerVector& groups_, const Rcpp::IntegerVector& unique_groups_, const Rcpp::NumericVector& group_weights_, const Rcpp::NumericVector& weights_, const Rcpp::NumericVector& weights_s_, const Rcpp::NumericVector& offset_, const Rcpp::NumericVector& offset_s_, const Rcpp::NumericVector& lambda_, const int nlambda, const double lambda_min_ratio, const double tau, const int maxit, const double tol, const int maxit_irls, const double tol_irls, const bool intercept_z, const bool intercept_s, const std::vector<std::string> penalty, const bool opposite_signs, const bool strongrule);
RcppExport SEXP _personalized2part_fit_twopart_cpp(SEXP X_SEXP, SEXP Z_SEXP, SEXP Xs_SEXP, SEXP S_SEXP, SEXP groups_SEXP, SEXP unique_groups_SEXP, SEXP group_weights_SEXP, SEXP weights_SEXP, SEXP weights_s_SEXP, SEXP offset_SEXP, SEXP offset_s_SEXP, SEXP lambda_SEXP, SEXP nlambdaSEXP, SEXP lambda_min_ratioSEXP, SEXP tauSEXP, SEXP maxitSEXP, SEXP tolSEXP, SEXP maxit_irlsSEXP, SEXP tol_irlsSEXP, SEXP intercept_zSEXP, SEXP intercept_sSEXP, SEXP penaltySEXP, SEXP opposite_signsSEXP, SEXP strongruleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type Z_(Z_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Xs_(Xs_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type S_(S_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type groups_(groups_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type unique_groups_(unique_groups_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type group_weights_(group_weights_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weights_(weights_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type weights_s_(weights_s_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type offset_(offset_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type offset_s_(offset_s_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< const int >::type nlambda(nlambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_min_ratio(lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit_irls(maxit_irlsSEXP);
    Rcpp::traits::input_parameter< const double >::type tol_irls(tol_irlsSEXP);
    Rcpp::traits::input_parameter< const bool >::type intercept_z(intercept_zSEXP);
    Rcpp::traits::input_parameter< const bool >::type intercept_s(intercept_sSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string> >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< const bool >::type opposite_signs(opposite_signsSEXP);
    Rcpp::traits::input_parameter< const bool >::type strongrule(strongruleSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_twopart_cpp(X_, Z_, Xs_, S_, groups_, unique_groups_, group_weights_, weights_, weights_s_, offset_, offset_s_, lambda_, nlambda, lambda_min_ratio, tau, maxit, tol, maxit_irls, tol_irls, intercept_z, intercept_s, penalty, opposite_signs, strongrule));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_personalized2part_fit_twopart_cpp", (DL_FUNC) &_personalized2part_fit_twopart_cpp, 24},
    {NULL, NULL, 0}
};

RcppExport void R_init_personalized2part(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
