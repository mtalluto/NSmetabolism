// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pressureCorrection
Rcpp::NumericVector pressureCorrection(Rcpp::NumericVector P, Rcpp::NumericVector elev, double newElev);
RcppExport SEXP _NSmetabolism_pressureCorrection(SEXP PSEXP, SEXP elevSEXP, SEXP newElevSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type P(PSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type elev(elevSEXP);
    Rcpp::traits::input_parameter< double >::type newElev(newElevSEXP);
    rcpp_result_gen = Rcpp::wrap(pressureCorrection(P, elev, newElev));
    return rcpp_result_gen;
END_RCPP
}
// idw_matrix
NumericVector idw_matrix(const NumericMatrix& vals, NumericVector dist, double pow);
RcppExport SEXP _NSmetabolism_idw_matrix(SEXP valsSEXP, SEXP distSEXP, SEXP powSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type pow(powSEXP);
    rcpp_result_gen = Rcpp::wrap(idw_matrix(vals, dist, pow));
    return rcpp_result_gen;
END_RCPP
}
// idw_river
NumericVector idw_river(List vals, List dist, List nbQ, NumericVector Q);
RcppExport SEXP _NSmetabolism_idw_river(SEXP valsSEXP, SEXP distSEXP, SEXP nbQSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type vals(valsSEXP);
    Rcpp::traits::input_parameter< List >::type dist(distSEXP);
    Rcpp::traits::input_parameter< List >::type nbQ(nbQSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(idw_river(vals, dist, nbQ, Q));
    return rcpp_result_gen;
END_RCPP
}
// computeInputDOFlux
double computeInputDOFlux(const Rcpp::NumericVector& upQ, const Rcpp::NumericVector& upDO, double latQ, double latDO);
RcppExport SEXP _NSmetabolism_computeInputDOFlux(SEXP upQSEXP, SEXP upDOSEXP, SEXP latQSEXP, SEXP latDOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type upQ(upQSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type upDO(upDOSEXP);
    Rcpp::traits::input_parameter< double >::type latQ(latQSEXP);
    Rcpp::traits::input_parameter< double >::type latDO(latDOSEXP);
    rcpp_result_gen = Rcpp::wrap(computeInputDOFlux(upQ, upDO, latQ, latDO));
    return rcpp_result_gen;
END_RCPP
}
// computeAdvection
double computeAdvection(double inputDOFlux, double DOconc, const Rcpp::NumericVector& data);
RcppExport SEXP _NSmetabolism_computeAdvection(SEXP inputDOFluxSEXP, SEXP DOconcSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type inputDOFlux(inputDOFluxSEXP);
    Rcpp::traits::input_parameter< double >::type DOconc(DOconcSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(computeAdvection(inputDOFlux, DOconc, data));
    return rcpp_result_gen;
END_RCPP
}
// computeGPP
Rcpp::NumericVector computeGPP(const Rcpp::NumericVector& PAR, const Rcpp::NumericVector& lP1, const Rcpp::NumericVector& lP2);
RcppExport SEXP _NSmetabolism_computeGPP(SEXP PARSEXP, SEXP lP1SEXP, SEXP lP2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type PAR(PARSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type lP1(lP1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type lP2(lP2SEXP);
    rcpp_result_gen = Rcpp::wrap(computeGPP(PAR, lP1, lP2));
    return rcpp_result_gen;
END_RCPP
}
// inSituER
Rcpp::NumericVector inSituER(const Rcpp::NumericVector& temperature, const Rcpp::NumericVector& ER24_20);
RcppExport SEXP _NSmetabolism_inSituER(SEXP temperatureSEXP, SEXP ER24_20SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type temperature(temperatureSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ER24_20(ER24_20SEXP);
    rcpp_result_gen = Rcpp::wrap(inSituER(temperature, ER24_20));
    return rcpp_result_gen;
END_RCPP
}
// computeRF
Rcpp::NumericVector computeRF(const Rcpp::NumericVector& temp, const Rcpp::NumericVector& pressure, const Rcpp::NumericVector& DO, const Rcpp::NumericVector& k600);
RcppExport SEXP _NSmetabolism_computeRF(SEXP tempSEXP, SEXP pressureSEXP, SEXP DOSEXP, SEXP k600SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type pressure(pressureSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type DO(DOSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type k600(k600SEXP);
    rcpp_result_gen = Rcpp::wrap(computeRF(temp, pressure, DO, k600));
    return rcpp_result_gen;
END_RCPP
}
// kT
Rcpp::NumericVector kT(const Rcpp::NumericVector& temp, const Rcpp::NumericVector& k600);
RcppExport SEXP _NSmetabolism_kT(SEXP tempSEXP, SEXP k600SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type k600(k600SEXP);
    rcpp_result_gen = Rcpp::wrap(kT(temp, k600));
    return rcpp_result_gen;
END_RCPP
}
// osat
Rcpp::NumericVector osat(const Rcpp::NumericVector& temp, const Rcpp::NumericVector& P);
RcppExport SEXP _NSmetabolism_osat(SEXP tempSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(osat(temp, P));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NSmetabolism_pressureCorrection", (DL_FUNC) &_NSmetabolism_pressureCorrection, 3},
    {"_NSmetabolism_idw_matrix", (DL_FUNC) &_NSmetabolism_idw_matrix, 3},
    {"_NSmetabolism_idw_river", (DL_FUNC) &_NSmetabolism_idw_river, 4},
    {"_NSmetabolism_computeInputDOFlux", (DL_FUNC) &_NSmetabolism_computeInputDOFlux, 4},
    {"_NSmetabolism_computeAdvection", (DL_FUNC) &_NSmetabolism_computeAdvection, 3},
    {"_NSmetabolism_computeGPP", (DL_FUNC) &_NSmetabolism_computeGPP, 3},
    {"_NSmetabolism_inSituER", (DL_FUNC) &_NSmetabolism_inSituER, 2},
    {"_NSmetabolism_computeRF", (DL_FUNC) &_NSmetabolism_computeRF, 4},
    {"_NSmetabolism_kT", (DL_FUNC) &_NSmetabolism_kT, 2},
    {"_NSmetabolism_osat", (DL_FUNC) &_NSmetabolism_osat, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_NSmetabolism(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
