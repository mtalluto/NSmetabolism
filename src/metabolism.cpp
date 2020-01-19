#include <Rcpp.h>
#include <cmath>
#include "../inst/include/metabolism.h"


long double NSM::log_sum_exp(long double v1, long double v2) {
	long double mv = std::max(v1, v2);
	return 	mv + std::log(std::exp(v1 - mv) + std::exp(v2 - mv));
}


double NSM::gpp(double PAR, double lP1, double lP2) {
	if(Rcpp::NumericVector::is_na(lP2))
		return NSM::gpp(PAR, lP1);

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

double NSM::er(double temp, double ER24_20) {
	if(ER24_20 > 0)
		throw std::range_error("DO concentration must be <= 0");
	return ER24_20 * pow(1.045, temp - 20);
}



double NSM::reaeration(double temp, double pressure, double DO, double k600) {
	// kT in m/day, osat in g/m^3, yields g/(m^2*day)
	return kT(temp, k600) * (osat(temp, pressure) - DO);
}


double NSM::kT(double temp, double k600) {
	if(k600 < 0)
		throw std::range_error("k600 must be positive");

	// compute Schmidt number for oxygen (dimensionless)
	// parameters from Wanninkhof 1992. appendix
	double Sc = 1800.6 - 120.10 * temp + 3.7818 * std::pow(temp,2) - 0.047608 * std::pow(temp, 3);
	
	// Van de Bogert et al eqn 5
	return k600 * std::pow(Sc / 600, -0.5);
}



double NSM::osat(double temp, double P) {  
	double hPaPerAtm = 1013.2500;

	double tempK = temp + 273.15;
	double Patm = P / hPaPerAtm;

	// C*o, the unit standard atmospheric concentration by volume for oxygen; in mg/L
	// eqn 32 from Benson and Krause 1984.
	double Cstaro = std::exp(-139.34411 + (1.575701e5 / tempK) - 
		(6.642308e7 / std::pow(tempK, 2)) + (1.243800e10 / std::pow(tempK, 3)) - 
		(8.621949e11 / std::pow(tempK, 4)));

	// Benson and Krause 1980 eqn 13, 
	// the negative of the second pressure coefficient in the virial expansion for 
	// the real behavior of oxygen.
	double theta = 0.000975 - 1.426e-5 * temp + 6.436e-8 * std::pow(temp, 2);

	// Benson and Krause 1980 eqn 23
	// the saturated vapor pressure of water in atmospheres at the temperature of equilibrium.
	// in atmospheres
	double Pwv = std::exp(11.8571 - (3840.7 / tempK) - (216961 / std::pow(tempK, 2)));

	// Benson and Krause 1980 eqn 28
	return Cstaro * ((Patm - Pwv) * (1 - theta * Patm)) / ((1 - Pwv) * (1 - theta));
}

