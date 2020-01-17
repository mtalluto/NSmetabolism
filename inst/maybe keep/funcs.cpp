#include <algorithm>
#include <cmath>
#include <Rcpp.h>
#include <vector>
#include "../inst/include/funcs.h"



//' @export
// [[Rcpp::export]]
double pressureCorrection (double P, double elev, double newElev) {
	double a = 2.25577e-5;
	double b = 5.25588;
	double PPa = P * 100; // convert from hPa
	double seaLevelP = PPa / std::pow(1 - a * elev, b);
	return((seaLevelP * std::pow(1 - a * newElev, b))/100);
}

