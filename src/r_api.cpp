#include <Rcpp.h>
#include <vector>
#include "../inst/include/r_api.h"
#include "../inst/include/check.h"
#include "../inst/include/do.h"

/*
	DISSOLVED OXYGEN FUNCTIONS
	these functions are contained in do.h
*/

//' Computes the input dissolved oxygen flux
//' @param upQ Upstream discharge (m^3/s)
//' @param upDO Upstream dissolved oxygen concentration (g/m^3)
//' @param latQ Lateral input discharge (m^3/s)
//' @param latDO Lateral input DO concentration (g/m^3)
//' @return double; the input flux (mg/min) of dissolved oxygen
//' 
// [[Rcpp::export]]
double computeInputDOFlux(const Rcpp::NumericVector &upQ, const Rcpp::NumericVector &upDO, 
	double latQ, double latDO) {
	if(upQ.size() != upDO.size())
		throw std::length_error("upQ and upDO must have the same length");

	std::vector<double> Q = Rcpp::as<std::vector<double> >(upQ);
	std::vector<double> DO = Rcpp::as<std::vector<double> >(upDO);
	Q.push_back(latQ);
	DO.push_back(latDO);

	check_DO(Rcpp::wrap(DO));
	check_Q(Rcpp::wrap(Q));

	return NSM::oxyMassFlux(DO, Q);
}

//' Compute advective flux in g/(m^3 * min)
//' @param inputDOFlux dissolved oxygen mass flux of input water (g/[m^3 * min])
//' @param DOconc dissolved oxygen concentration (g/m^3)
//' @param data Fixed site data, including at least Q, area, and dx
//' @return dissolved oxygen concentration flux (g/[m^3 * min])
//' 
// [[Rcpp::export]]
double computeAdvection(double inputDOFlux, double DOconc, const Rcpp::NumericVector &data) {
	check_DO(Rcpp::NumericVector {inputDOFlux, DOconc});
	check_site_data(data);

	return NSM::advection(inputDOFlux, DOconc, data["Q"], data["area"], data["dx"]);
}


//' Compute gross primary productivity from light (g O2 / [m^2 * day])
//' @details P1 is in units of (W*day)/gO2, P2 is in (m^2 * day) / gO2
//' @param lP1 The log of the slope of the PI curve
//' @param lP2 The log of the saturation term of the PI curve
//' @param data Fixed site data, including at least Q, area, and dx
//' @references Uehlinger U. et al. 2000. Variability of photosynthesis-irradiance curves and
//' 		ecosystem respiration in a small river. Freshwater Biology 44:493-507.
//' @return double; computed GPP 
//' 
// [[Rcpp::export]]
Rcpp::NumericVector computeGPP(const Rcpp::NumericVector &PAR, const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2) {
	check_light(PAR);
	if(PAR.size() != lP1.size() || PAR.size() != lP2.size())
		throw std::length_error("All GPP inputs must have the same length");

	Rcpp::NumericVector result (PAR.size());
	for(int i = 0; i < PAR.size(); ++i) {
		if(Rcpp::NumericVector::is_na(lP2(i)))
			result(i) = NSM::gpp(PAR(i), lP1(i));
		else
			result(i) = NSM::gpp(PAR(i), lP1(i), lP2(i));
	}
	return result;
}


//' Compute ecosystem respiration at in-situ temperature (gO2 / [m^2 * day])
//' @param temperature Water temperature (degrees C)
//' @param ER24_20 Ecosystem respiration rate at 20 degrees (gO2 / [m^2 * day])
//' @return Temperature-adjusted ER
//' 
// [[Rcpp::export]]
Rcpp::NumericVector inSituER(const Rcpp::NumericVector &temperature, 
			const Rcpp::NumericVector &ER24_20) {

	if(ER24_20.size() != 1 && ER24_20.size() != temperature.size())
		throw std::length_error("ER24_20 must have length equal to 1 or to length(temperature)");


	Rcpp::NumericVector result (temperature.size());

	double er24 = ER24_20(0);
	for(int i = 0; i < temperature.size(); ++i) {
		if(ER24_20.size() > 1)
			er24 = ER24_20(i);
		result(i) = NSM::er(temperature(i), er24);
	}
	return result;
}

//' Computes reaeration flux (g/(m^2 * day))
//' 
//' @param temp Water temperature (degrees C)
//' @param pressure Atmospheric pressure (hPa)
//' @param DO dissolved oxygen concentration, mg/L
//' @param k600 Gas transfer coefficient for Schmidt number of 600
//' @return computed rearation flux
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector computeRF(const Rcpp::NumericVector &temp, 
		const Rcpp::NumericVector &pressure, const Rcpp::NumericVector &DO, 
		const Rcpp::NumericVector &k600) {
	// check vector dimensions
	if(temp.size() != pressure.size() || temp.size() != DO.size())
		throw std::length_error("temp, pressure, and DO must all be the same length");
	if(temp.size() != k600.size() && k600.size() != 1)
		throw std::length_error("Length of k600 must be 1 or the same as other vectors");
	check_pressure(pressure);
	check_DO(DO);
	check_k(k600);

	Rcpp::NumericVector result (temp.size());
	double k = k600(0);
	for(int i = 0; i < temp.size(); ++i) {
		if(k600.size() > 1)
			k = k600(i);
		result(i) = NSM::reaeration(temp(i), pressure(i), DO(i), k);
	}
	return result;
}


//' Compute gas transfer velocity for oxygen at a given temperature (m/day)
//' @param temp Water temperature (degrees C)
//' @param k600 Gas transfer coefficient for Schmidt number of 600 (m/day)
//' @references Wanninkhof R. (1992). Relationship between wind speed and gas exchange over the
//'    ocean. Journal of Geophysical Research, 97, 7373.\n
//'    Van de Bogert, M.C., Carpenter, S.R., Cole, J.J. & Pace, M.L. (2007). Assessing pelagic 
//'    and benthic metabolism using free water measurements. Limnology and Oceanography: 
//'    Methods, 5, 145–155.
//' @return k at the given temperature (m/day)
// [[Rcpp::export]]
Rcpp::NumericVector kT(const Rcpp::NumericVector &temp, const Rcpp::NumericVector &k600) {
	if(temp.size() != k600.size())
		throw std::length_error("kT input vectors must have the same length");

	Rcpp::NumericVector result (temp.size());
	for(int i = 0; i < temp.size(); ++i) {
		result(i) = NSM::kT(temp(i), k600(i));
	}
	return result;
}



//' Compute dissolved oxygen saturation given temperature and pressure (mg/L OR g/m^3)
//' @param temp Water temperature (degrees C)
//' @param P atmospheric pressure, in hPa
//' @references Benson BB and Krause D. 1984. The concentration and isotopic fractionation of
//'     oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere. 
//'     Limnol. Oceanogr., 29, 620–632.
//'     
//' Benson BB and Krause D. 1980. The concentration and isotopic fractionation of gases
//'     		dissolved in freshwater in equilibrium with the atmosphere. 1. Oxygen.  Limnol. 
//'     		Oceanogr., 25(4):662-671
//' @return oxygen saturation concentration at given temperature and pressure
// [[Rcpp::export]]
Rcpp::NumericVector osat(const Rcpp::NumericVector &temp, const Rcpp::NumericVector &P) {
	if(temp.size() != P.size())
		throw std::length_error("osat input vectors must have the same length");
	check_pressure(P);

	Rcpp::NumericVector result (temp.size());
	for(int i = 0; i < temp.size(); ++i) {
		result(i) = NSM::osat(temp(i), P(i));
	}
	return result;
}
