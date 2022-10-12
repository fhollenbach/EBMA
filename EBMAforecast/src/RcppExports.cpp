// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// oneMultinomCalt
IntegerVector oneMultinomCalt(NumericVector probs);
RcppExport SEXP _EBMAforecast_oneMultinomCalt(SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(oneMultinomCalt(probs));
    return rcpp_result_gen;
END_RCPP
}
// getRGamma
NumericVector getRGamma(double shape);
RcppExport SEXP _EBMAforecast_getRGamma(SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(getRGamma(shape));
    return rcpp_result_gen;
END_RCPP
}
// isNA
LogicalVector isNA(NumericVector x);
RcppExport SEXP _EBMAforecast_isNA(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(isNA(x));
    return rcpp_result_gen;
END_RCPP
}
// emNorm
Rcpp::List emNorm(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericMatrix RSQ, Rcpp::NumericVector W, double tol, int maxIter, double wisdom, double sigma2);
RcppExport SEXP _EBMAforecast_emNorm(SEXP outcomeSEXP, SEXP predictionSEXP, SEXP RSQSEXP, SEXP WSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP wisdomSEXP, SEXP sigma2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type RSQ(RSQSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type wisdom(wisdomSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    rcpp_result_gen = Rcpp::wrap(emNorm(outcome, prediction, RSQ, W, tol, maxIter, wisdom, sigma2));
    return rcpp_result_gen;
END_RCPP
}
// GibbsNormal
Rcpp::List GibbsNormal(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, Rcpp::NumericVector alpha, double sigma, int iterations, int burnin, int thin);
RcppExport SEXP _EBMAforecast_GibbsNormal(SEXP outcomeSEXP, SEXP predictionSEXP, SEXP WSEXP, SEXP alphaSEXP, SEXP sigmaSEXP, SEXP iterationsSEXP, SEXP burninSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type W(WSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(GibbsNormal(outcome, prediction, W, alpha, sigma, iterations, burnin, thin));
    return rcpp_result_gen;
END_RCPP
}
// GibbsNormalMissing
Rcpp::List GibbsNormalMissing(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, Rcpp::NumericVector alpha, double sigma, int iterations, int burnin, int thin);
RcppExport SEXP _EBMAforecast_GibbsNormalMissing(SEXP outcomeSEXP, SEXP predictionSEXP, SEXP WSEXP, SEXP alphaSEXP, SEXP sigmaSEXP, SEXP iterationsSEXP, SEXP burninSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type W(WSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(GibbsNormalMissing(outcome, prediction, W, alpha, sigma, iterations, burnin, thin));
    return rcpp_result_gen;
END_RCPP
}
// emLogit
Rcpp::List emLogit(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, double tol, int maxIter, double wisdom);
RcppExport SEXP _EBMAforecast_emLogit(SEXP outcomeSEXP, SEXP predictionSEXP, SEXP WSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP wisdomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type wisdom(wisdomSEXP);
    rcpp_result_gen = Rcpp::wrap(emLogit(outcome, prediction, W, tol, maxIter, wisdom));
    return rcpp_result_gen;
END_RCPP
}
// GibbsLogit
Rcpp::List GibbsLogit(Rcpp::NumericVector outcome, Rcpp::NumericMatrix prediction, Rcpp::NumericVector W, int iterations, int burnin, int thin);
RcppExport SEXP _EBMAforecast_GibbsLogit(SEXP outcomeSEXP, SEXP predictionSEXP, SEXP WSEXP, SEXP iterationsSEXP, SEXP burninSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type W(WSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(GibbsLogit(outcome, prediction, W, iterations, burnin, thin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EBMAforecast_oneMultinomCalt", (DL_FUNC) &_EBMAforecast_oneMultinomCalt, 1},
    {"_EBMAforecast_getRGamma", (DL_FUNC) &_EBMAforecast_getRGamma, 1},
    {"_EBMAforecast_isNA", (DL_FUNC) &_EBMAforecast_isNA, 1},
    {"_EBMAforecast_emNorm", (DL_FUNC) &_EBMAforecast_emNorm, 8},
    {"_EBMAforecast_GibbsNormal", (DL_FUNC) &_EBMAforecast_GibbsNormal, 8},
    {"_EBMAforecast_GibbsNormalMissing", (DL_FUNC) &_EBMAforecast_GibbsNormalMissing, 8},
    {"_EBMAforecast_emLogit", (DL_FUNC) &_EBMAforecast_emLogit, 6},
    {"_EBMAforecast_GibbsLogit", (DL_FUNC) &_EBMAforecast_GibbsLogit, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_EBMAforecast(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
