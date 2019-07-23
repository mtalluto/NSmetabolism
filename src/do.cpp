/*
	functions needed to compute dissolved oxygen concentration
*/

#include <Rcpp.h>
#include <cmath>
#include "../inst/include/funcs.h"

double dDOdt (const Rcpp::NumericVector &params, const Rcpp::NumericVector &data, 
			double inputDOMass, double DOPrev, double light, double waterTemp, double pressure);
double computeER(double temp, double ER24_20);
double computeGPP(double PAR, double lP1, double lP2);
double computeGPP_linear(double PAR, double lP1);
double computeRF(double temp, double pressure, double DO, double k600);
double kT(double temp, double k600);
double osat(double temp, double P);
double computeAdvection(double inputDOMass, double outputDO, double Q, double area, double dx);


// idea - input checking is all over the place
// only do input checking on r entry points
// so need an r_api file that has r entry points to various functions
// these do input checking and then call c++ functions
// then I can have a check_params() and check_data() function that handles this

/**
	* Compute derivative of dissolved oxygen with respect to time
	*
	* @details Expected named parameters in the model are:
		* `k600`: the gas exchange coefficient
		* `lP1`: log of the slope of the PI curve
		* `lP2`: log of the saturation term of the PI curve
		* `ER24`: ER for a 24 h period at 20 degrees

	Site specific data necessary for dDOdt are:
		* `Q`: discharge (m^3 / sec)
		* `area`: Cross-sectional area (m^2)
		* `dx` length of stream segment (meters)
		* `z` stream depth (meters)

	* @param params Named vector of model parameters
	* @param data Named vector of site- (i.e., pixel-) specific data
	* @param inputDOMass DO mass from upstream and lateral sources (mg)
	* @param DOPrev DO concentration at previuos time step (mg/L)
	* @param light Light irradiance (W/m^2)
	* @param waterTemp Water temperature (degrees C)
	* @param pressure Atmospheric pressure (hPa)
	* @return double; derivative of dissolved oxygen with respect to time
*/
// [[Rcpp::export]]
double dDOdt (const Rcpp::NumericVector &params, const Rcpp::NumericVector &data, 
			double inputDOMass, double DOPrev, double light, double waterTemp, double pressure) {
	if(data['z'] <= 0)
		throw std::range_error("Depth (z) must be positive");

	double advection = computeAdvection(inputDOMass, DOPrev, data["Q"], data["area"], data["dx"]);
	double gpp = computeGPP(light, params["lP1"], params["lP2"]);
	double er = computeER(waterTemp, params["ER24"]);
	double rf = computeRF(waterTemp, pressure, DOPrev, params["k600"]);

	double ddodt = advection + (gpp + er + rf) / data["z"];
	return ddodt;
}


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
	if(PAR < 0)
		throw std::range_error("Light must be >= 0");
	if(PAR == 0)
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
double computeRF(double temp, double pressure, double DO, double k600) {
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
	* Compute dissolved oxygen saturation given temperature and pressure
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


/**
	* Compute transport component
	* @param inputDOMass dissolved oxygen mass of input water
	* @param outputDO dissolved oxygen concentration of focal pixel
	* @param Q discharge
	* @param area cross sectional area
	* @param dx length of pixel
	* @return transported DO mass
*/
// [[Rcpp::export]]
double computeAdvection(double inputDOMass, double outputDO, double Q, double area, double dx) {
	if(inputDOMass < 0 || outputDO < 0)
		throw std::range_error("inputDOMass and outputDO must be >= 0");
	if(Q < 0)
		throw std::range_error("Q must be >= 0");
	if(area <= 0 || dx <= 0)
		throw std::range_error("Area and dx must positive");

	double outputDOMass = Q * outputDO;
	return (-1/area) * (outputDOMass - inputDOMass)/dx;
}
