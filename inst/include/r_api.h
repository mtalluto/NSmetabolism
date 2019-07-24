#ifndef RAPI_H
#define RAPI_H
#include <Rcpp.h>

double computeInputDOFlux(const Rcpp::NumericVector &upQ, const Rcpp::NumericVector &upDO, 
	double latQ, double latDO);

double computeAdvection(double inputDOFlux, double DOconc, const Rcpp::NumericVector &data);

Rcpp::NumericVector computeGPP(const Rcpp::NumericVector &PAR, const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2);

Rcpp::NumericVector inSituER(const Rcpp::NumericVector &temperature, 
			const Rcpp::NumericVector &ER24_20);

Rcpp::NumericVector computeRF(const Rcpp::NumericVector &temp, 
		const Rcpp::NumericVector &pressure, const Rcpp::NumericVector &DO, 
		const Rcpp::NumericVector &k600);

Rcpp::NumericVector kT(const Rcpp::NumericVector &temp, const Rcpp::NumericVector &k600);

Rcpp::NumericVector osat(const Rcpp::NumericVector &temp, const Rcpp::NumericVector &P);

#endif