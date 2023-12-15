// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// glarma_cpp
List glarma_cpp(const arma::mat& s, double rho, double mtol, int maxIterations, double ltol);
RcppExport SEXP _Glarmadillo_glarma_cpp(SEXP sSEXP, SEXP rhoSEXP, SEXP mtolSEXP, SEXP maxIterationsSEXP, SEXP ltolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type mtol(mtolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIterations(maxIterationsSEXP);
    Rcpp::traits::input_parameter< double >::type ltol(ltolSEXP);
    rcpp_result_gen = Rcpp::wrap(glarma_cpp(s, rho, mtol, maxIterations, ltol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Glarmadillo_glarma_cpp", (DL_FUNC) &_Glarmadillo_glarma_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_Glarmadillo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
