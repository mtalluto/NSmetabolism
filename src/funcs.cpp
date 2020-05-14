#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include <vector>
#include "../inst/include/funcs.h"

long double log_sum_exp(long double v1, long double v2) {
	long double mv = std::max(v1, v2);
	return 	mv + std::log(std::exp(v1 - mv) + std::exp(v2 - mv));
}

//' Compute pressure at a given elevation from a calibration point
//' @param P vector of pressure observations (e.g., a time series)
//' @param elev Vector of elevations (in meters)
//' @param calibE Elevaiton of calibration point
//' @return Pressure in hPa; in the form of a matrix with elevation as rows and pressure as columns
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pressureCorrection (Rcpp::NumericVector P, Rcpp::NumericVector elev, 
		double newElev) {
	if(P.size() != elev.size())
		Rcpp::stop("P and elev must be the same length");
	Rcpp::NumericVector result;
	for(int i = 0; i < P.size(); ++i)
		result.push_back(pressureCorrection(P[i], elev[i], newElev));
	return result;
}

double pressureCorrection (double P, double elev, double newElev) {
	double a = 2.25577e-5;
	double b = 5.25588;
	double PPa = P * 100; // convert from hPa
	double seaLevelP = PPa / std::pow(1 - a * elev, b);
	return((seaLevelP * std::pow(1 - a * newElev, b))/100);
}