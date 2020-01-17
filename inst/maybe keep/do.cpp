/*
	functions needed to compute dissolved oxygen concentration
*/

#include "../inst/include/funcs.h"

// double dDOdt (const Rcpp::NumericVector &params, const Rcpp::NumericVector &data, 
// 			double inputDOMass, double DOPrev, double light, double waterTemp, double pressure);
// double computeER(double temp, double ER24_20);
// double computeGPP(double PAR, double lP1, double lP2);
// double computeGPP_linear(double PAR, double lP1);
// double computeRF(double temp, double pressure, double DO, double k600);
// double kT(double temp, double k600);
// double osat(double temp, double P);
// double computeAdvection(double inputDOMass, double outputDO, double Q, double area, double dx);



double NSM::oxyMassFlux (const std::vector<double> &DO, const std::vector<double> &Q) {
	double result = 0;
	for(int i = 0; i < Q.size(); ++i)
		result += oxyMassFlux(DO.at(i), Q.at(i));
	return result;
}


double NSM::oxyMassFlux (double DO, double Q) {
	// converts discharge into m^3/min
	double sPerMin = 60;
	Q *= sPerMin;
	return Q * DO;
}


double NSM::advection(double inputFlux, double DOconc, double Q, double area, double dx) {
	double outputFlux = NSM::oxyMassFlux(DOconc, Q);
	return (-1/area) * (outputFlux - inputFlux)/dx;
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
// double dDOdt (const Rcpp::NumericVector &params, const Rcpp::NumericVector &data, 
// 			double inputDOMass, double DOPrev, double light, double waterTemp, double pressure) {
// 	if(data['z'] <= 0)
// 		throw std::range_error("Depth (z) must be positive");

// 	double advection = computeAdvection(inputDOMass, DOPrev, data["Q"], data["area"], data["dx"]);
// 	NOTE UNIT CHANGE - advection in g/(m^3 * min), others in g/(m^3 * day)
// 	double gpp = computeGPP(light, params["lP1"], params["lP2"]);
// 	double er = computeER(waterTemp, params["ER24"]);
// 	double rf = computeRF(waterTemp, pressure, DOPrev, params["k600"]);

// convert everything to min
// 	double ddodt = advection + ((gpp + er + rf) / data["z"]) / (24 * 60); 
// 	return ddodt;
// }



















