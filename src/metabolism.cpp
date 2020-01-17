#include <Rcpp.h>
#include <cmath>
#include "../inst/include/metabolism.h"


long double NSM::log_sum_exp(long double v1, long double v2) {
	long double mv = std::max(v1, v2);
	return 	mv + std::log(std::exp(v1 - mv) + std::exp(v2 - mv));
}


double NSM::gpp(double PAR, double lP1, double lP2) {
	if(PAR == 0)
		return 0;

	// Uehlinger et al 2000 eq 3b
	// GPP = PAR/(P1 + P2 * PAR)
	double lpar = std::log(PAR);
	double lGPP = lpar - NSM::log_sum_exp(lP1, lP2 + lpar);
	return std::exp(lGPP);
}


double NSM::gpp(double PAR, double lP1) {
	if(PAR == 0)
		return 0;

	double lGPP = std::log(PAR) - lP1;
	return std::exp(lGPP);
}
