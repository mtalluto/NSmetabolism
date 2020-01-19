#ifndef RAPI_H
#define RAPI_H
#include <Rcpp.h>

double computeInputDOFlux(const Rcpp::NumericVector &upQ, const Rcpp::NumericVector &upDO, 
	double latQ, double latDO);

double computeAdvection(double inputDOFlux, double DOconc, const Rcpp::NumericVector &data);





#endif