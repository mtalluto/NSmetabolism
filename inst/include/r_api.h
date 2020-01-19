#ifndef rapi_h
#define rapi_h
#include <Rcpp.h>
#include <string>

// compute GPP given a a vector of light values and parameters
Rcpp::NumericVector computeGPP(const Rcpp::NumericVector &PAR, const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2);

Rcpp::NumericVector inSituER(const Rcpp::NumericVector &temperature, 
			const Rcpp::NumericVector &ER24_20);

Rcpp::NumericVector computeRF(const Rcpp::NumericVector &temp, 
		const Rcpp::NumericVector &pressure, const Rcpp::NumericVector &DO, 
		const Rcpp::NumericVector &k600);

Rcpp::NumericVector kT(const Rcpp::NumericVector &temp, const Rcpp::NumericVector &k600);

Rcpp::NumericVector osat(const Rcpp::NumericVector &temp, const Rcpp::NumericVector &P);



/*
	Pixel-based entry points
*/

// compute metabolism from pixel data
Rcpp::NumericMatrix pixelMetabolism(const Rcpp::DataFrame &pixdf,  
			const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature, 
			const Rcpp::NumericMatrix &pressure, const Rcpp::NumericVector &do_init,
			const Rcpp::NumericVector &lP1, const Rcpp::NumericVector &lP2, 
			const Rcpp::NumericVector &er24_20, const Rcpp::NumericVector &k600,
			std::string variable);


#endif