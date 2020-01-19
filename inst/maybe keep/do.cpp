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
double dDOdt (const Rcpp::NumericVector &params, const Rcpp::NumericVector &data, 
			double inputDOMass, double DOPrev, double light, double waterTemp, double pressure) {
	if(data['z'] <= 0)
		throw std::range_error("Depth (z) must be positive");

	double advection = computeAdvection(inputDOMass, DOPrev, data["Q"], data["area"], data["dx"]);
	NOTE UNIT CHANGE - advection in g/(m^3 * min), others in g/(m^3 * day)
	double gpp = computeGPP(light, params["lP1"], params["lP2"]);
	double er = computeER(waterTemp, params["ER24"]);
	double rf = computeRF(waterTemp, pressure, DOPrev, params["k600"]);

convert everything to min
	double ddodt = advection + ((gpp + er + rf) / data["z"]) / (24 * 60); 
	return ddodt;
}



















