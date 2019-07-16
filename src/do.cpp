/*
	functions needed to compute dissolved oxygen concentration
*/

#include <Rcpp.h>
#include <cmath>
#include "../inst/include/funcs.h"

double computeGPP(double PAR, double lP1, double lP2);
double computeGPP_linear(double PAR, double lP1);
double compute_rf(double temp, double pressure, double DO, double k600);
double kT(double temp, double k600);
double osat(double temp, double P);


/**
  * Compute ecosystem respiration at a given temperature
  * @param temp Water temperature (degrees C)
  * @param ER24_20 Ecosystem respiration rate for a 24-hour period at 20 degrees
*/
// [[Rcpp::export]]
double computeER(double temp, double ER24_20) {
	if(ER24_20 > 0)
		throw std::range_error("DO concentration must be <= 0");
	return (ER24_20 / (60*24)) * pow(1.045, temp - 20);
}

/**
  * Compute GPP from light
  * @param PAR Photosynthetically active radiation, in W/m^2 
  * @param lP1 The log of the slope of the PI curve
  * @param lP2 The log of the saturation term of the PI curve
*/
// [[Rcpp::export]]
double computeGPP(double PAR, double lP1, double lP2) {
	if(PAR <= 0)
		return 0;

	// Uehlinger et al 2000 eq 3b
	// GPP = PAR/(P1 + P2 * PAR)
	double lpar = std::log(PAR);
	double lGPP = lpar - log_sum_exp(lP1, lP2 + lpar);
	return std::exp(lGPP);
}


/**
  * Compute GPP from light
  * @param PAR Photosynthetically active radiation, in W/m^2 
  * @param lP1 The log of the slope of the PI curve
  * @param lP2 The log of the saturation term of the PI curve
*/
// [[Rcpp::export]]
double computeGPP_linear(double PAR, double lP1) {
	if(PAR <= 0)
		return 0;

	// Uehlinger et al 2000 eq 3b
	// GPP = PAR/(P1 + P2 * PAR)
	double lGPP = std::log(PAR) - lP1;
	return std::exp(lGPP);
}


/**
  * Computes reaeration flux
  *
  * @param temp Water temperature (degrees C)
  * @param pressure Atmospheric pressure (hPa)
  * @param DO dissolved oxygen concentration, mg/L
  * @param k600 Gas transfer coefficient for Schmidt number of 600
  * @return computed rearation flux
  * @export
*/
// [[Rcpp::export]]
double compute_rf(double temp, double pressure, double DO, double k600) {
	if(DO < 0)
		throw std::range_error("DO concentration must be positive");

	return kT(temp, k600) * (osat(temp, pressure) - DO);
}

/**
  * Compute gas transfer velocity for oxygen at a given temperature
  * @param temp Water temperature (degrees C)
  * @param k600 Gas transfer coefficient for Schmidt number of 600
  * @references Wanninkhof R. (1992). Relationship between wind speed and gas exchange over the
  *    ocean. Journal of Geophysical Research, 97, 7373.\n
  *    Van de Bogert, M.C., Carpenter, S.R., Cole, J.J. & Pace, M.L. (2007). Assessing pelagic 
  *    and benthic metabolism using free water measurements. Limnology and Oceanography: 
  *    Methods, 5, 145–155.
  * @return k at the given temperature
*/
double kT(double temp, double k600) {
	if(k600 < 0)
		throw std::range_error("k600 must be positive");

	// compute Schmidt number for oxygen
	// parameters from Wanninkhof 1992. appendix
	double Sc = 1800.6 - 120.10 * temp + 3.7818 * std::pow(temp,2) - 0.047608 * std::pow(temp, 3);
	
	// Van de Bogert et al eqn 5
	return k600 * std::pow(Sc / 600, -0.5);
}


/**
  * Compute oxygen saturation
  * @param temp Water temperature (degrees C)
  * @param P atmospheric pressure, in hPa
  * references Benson BB and Krause D. 1984. The concentration and isotopic fractionation of
  *     oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. 
  *     Limnol. Oceanogr., 29, 620–632.
  * @return oxygen saturation concentration at given temperature and pressure
*/
double osat(double temp, double P) {  
	if(P < 0)
		throw std::range_error("pressure must be positive");

	double hPaPerAtm = 1013.2500;

	double tempK = temp + 273.15;
	double Patm = P / hPaPerAtm;

	// C*o, the unit standard atmospheric concentration by volume for oxygen; in mg/kg
	// eqn 31 from Benson and Krause 1984.
	double Cstaro = std::exp(-1.3874202e2 + (1.572288e5 / tempK) - 
		(6.637149e7 / std::pow(tempK, 2)) + (1.243678e10 / std::pow(tempK, 3)) - 
		(8.621061e11 / std::pow(tempK, 4)));

	// eqn 13, 
	// the negative of the second pressure coefficient in the virial expansion for  gas 
	// the real behavior of oxygen.
	double theta = 0.000975 - 1.426e-5 * temp + 6.436e-8 * std::pow(temp, 2);

	// eqn 23
	// the saturated vapor pressure of water in atmospheres at the temperature of equilibrium.
	// in atmospheres
	double Pwv = std::exp(11.8571 - (3840.7 / tempK) - (216961 / std::pow(tempK, 2)));

	// eqn 28
	return Cstaro * ((Patm - Pwv) * (1 - theta * Patm)) / ((1 - Pwv) * (1 - theta));
}
