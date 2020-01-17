#ifndef rapi_h
#define rapi_h
#include <Rcpp.h>

// compute GPP given a a vector of light values and parameters
Rcpp::NumericVector computeGPP(const Rcpp::NumericVector &PAR, const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2);

// compute GPP from a data frame of pixel data
Rcpp::NumericMatrix pixelGPP(const Rcpp::DataFrame &pix, const Rcpp::NumericMatrix &light, 
		double lP1, double lP2);

#endif