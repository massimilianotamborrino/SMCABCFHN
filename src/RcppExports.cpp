// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// my_fun
NumericMatrix my_fun(int N, NumericVector mu, NumericMatrix Sigma);
RcppExport SEXP _SMCABCFHN_my_fun(SEXP NSEXP, SEXP muSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(my_fun(N, mu, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// FHN_prior2_
NumericVector FHN_prior2_(NumericVector theta, int draw);
RcppExport SEXP _SMCABCFHN_FHN_prior2_(SEXP thetaSEXP, SEXP drawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type draw(drawSEXP);
    rcpp_result_gen = Rcpp::wrap(FHN_prior2_(theta, draw));
    return rcpp_result_gen;
END_RCPP
}
// FHN_prior_
NumericVector FHN_prior_(NumericVector theta, int draw);
RcppExport SEXP _SMCABCFHN_FHN_prior_(SEXP thetaSEXP, SEXP drawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type draw(drawSEXP);
    rcpp_result_gen = Rcpp::wrap(FHN_prior_(theta, draw));
    return rcpp_result_gen;
END_RCPP
}
// FHN_model_
NumericMatrix FHN_model_(NumericVector theta, double delta, NumericVector X0, int N);
RcppExport SEXP _SMCABCFHN_FHN_model_(SEXP thetaSEXP, SEXP deltaSEXP, SEXP X0SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(FHN_model_(theta, delta, X0, N));
    return rcpp_result_gen;
END_RCPP
}
// FHN_model_check_
List FHN_model_check_(NumericVector theta, double delta, NumericVector X0, int N);
RcppExport SEXP _SMCABCFHN_FHN_model_check_(SEXP thetaSEXP, SEXP deltaSEXP, SEXP X0SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(FHN_model_check_(theta, delta, X0, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SMCABCFHN_my_fun", (DL_FUNC) &_SMCABCFHN_my_fun, 3},
    {"_SMCABCFHN_FHN_prior2_", (DL_FUNC) &_SMCABCFHN_FHN_prior2_, 2},
    {"_SMCABCFHN_FHN_prior_", (DL_FUNC) &_SMCABCFHN_FHN_prior_, 2},
    {"_SMCABCFHN_FHN_model_", (DL_FUNC) &_SMCABCFHN_FHN_model_, 4},
    {"_SMCABCFHN_FHN_model_check_", (DL_FUNC) &_SMCABCFHN_FHN_model_check_, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SMCABCFHN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
