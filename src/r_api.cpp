#include <Rcpp.h>
#include <vector>
#include <memory>
#include "../inst/include/r_api.h"
#include "../inst/include/metabolism.h"
#include "../inst/include/check.h"
#include "../inst/include/params.h"
#include "../inst/include/pixel.h"


//' Compute gross primary productivity from light (g O2 / [m^2 * day])
//' @details P1 is in units of (W*day)/gO2, P2 is in (m^2 * day) / gO2
//' @param lP1 The log of the slope of the PI curve
//' @param lP2 The log of the saturation term of the PI curve
//' @param data Fixed site data, including at least Q, area, and dx
//' @references Uehlinger U. et al. 2000. Variability of photosynthesis-irradiance curves and
//' 		ecosystem respiration in a small river. Freshwater Biology 44:493-507.
//' @return double; computed GPP 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector computeGPP(const Rcpp::NumericVector &light, const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2) {
	check_light(light);
	if(light.size() != lP1.size() || light.size() != lP2.size())
		throw std::length_error("All GPP inputs must have the same length");

	Rcpp::NumericVector result (light.size());
	for(int i = 0; i < light.size(); ++i) {
		if(Rcpp::NumericVector::is_na(lP2(i)))
			result(i) = NSM::gpp(light(i), lP1(i));
		else
			result(i) = NSM::gpp(light(i), lP1(i), lP2(i));
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
//'    ocean. Journal of Geophysical Research, 97, 7373.
//'
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






//' @name pixelMetabolism
//' @aliases pixelGPP
//' @aliases pixelER
//' @title Compute metabolic parameters using pixels
//' @param pixdf The pixel data frame to construct the pixels, see 'details'
//' @param light Matrix of light readings; pixels in rows, light readings in columns
//' @param temperature Water temperature (degrees C)
//' @param lP1 The log of the slope of the PI curve
//' @param lP2 The log of the saturation term of the PI curve
//' @param er24_20 Ecosystem respiration rate at 20 degrees (gO2 / [m^2 * day])
//'	@param variable; the variable to compute
//'	@param dt; the size of a time step, in minutes
//' @details These functions all take the same signature, and compute the time series of
//'		GPP and ER for pixelGPP and pixelER, or run the full simulation for pixelMetabolism.
//'
//' The pixel df must have the following named elements
//'		* depth: depth of the pixel, in m
//'		* DO_i: initial dissolved oxygen concentration
//'
//' @return A named List with elements DO (a matrix with predicted DO time series), 
//'		GPP (daily GPP), and ER (daily in-situ ER)
// [[Rcpp::export]]
Rcpp::List pixelMetabolism(const Rcpp::DataFrame &pixdf,  
			const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature, 
			const Rcpp::NumericMatrix &pressure, const Rcpp::NumericVector &lP1, 
			const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
			const Rcpp::NumericVector &k600, double dt = 10) {

	std::vector<NSM::param_ptr> par_vector = NSM::param_from_r(lP1, lP2, er24_20, k600);
	std::vector<NSM::Pixel> pixes = 
		NSM::dfToPixel(pixdf, light, temperature, pressure, par_vector, dt);

	// compute metab for all
	Rcpp::NumericMatrix dailyGPP (light.nrow(), pixes.at(0).ndays());
	Rcpp::NumericMatrix dailyER (light.nrow(), pixes.at(0).ndays());
	Rcpp::NumericMatrix dissOx (light.nrow(), light.ncol());
	for(int p = 0; p < pixes.size(); ++p) {
		NSM::Pixel &pix = pixes.at(p);
		pix.simulate_os();
		for(int i = 0; i < pix.do_history().size(); ++i)
			dissOx(p, i) = pix.do_history().at(i);

		for(int d = 0; d < pix.ndays(); ++d) {
			dailyGPP(p, d) = pix.daily_gpp().at(d);
			dailyER(p, d) = pix.daily_er().at(d);
		}
	}

	Rcpp::List result (Rcpp::List::create(Rcpp::Named("DO") = dissOx, 
		Rcpp::Named("GPP") = dailyGPP, Rcpp::Named("ER") = dailyER));
	return result;
}

//' @rdname pixelMetabolism
//' @return Matrix; computed GPP for each pixel
// [[Rcpp::export]]
Rcpp::NumericMatrix pixelGPP(const Rcpp::DataFrame &pixdf,  
			const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature, 
			const Rcpp::NumericMatrix &pressure, const Rcpp::NumericVector &lP1, 
			const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
			const Rcpp::NumericVector &k600, double dt = 10) {

	std::vector<NSM::param_ptr> par_vector = NSM::param_from_r(lP1, lP2, er24_20, k600);
	std::vector<NSM::Pixel> pixes = 
		NSM::dfToPixel(pixdf, light, temperature, pressure, par_vector, dt);

	Rcpp::NumericMatrix result (light.nrow(), light.ncol());
	for(int t = 0; t < light.ncol(); ++t) {
		for(int p = 0; p < pixes.size(); ++p) {
			result(p, t) = pixes.at(p).gpp(t);
		}
	}
	return result;
}

//' @rdname pixelMetabolism
//' @return Matrix; computed ER for each pixel
// [[Rcpp::export]]
Rcpp::NumericMatrix pixelER(const Rcpp::DataFrame &pixdf,  
			const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature, 
			const Rcpp::NumericMatrix &pressure, const Rcpp::NumericVector &lP1, 
			const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
			const Rcpp::NumericVector &k600, double dt = 10) {
	
	std::vector<NSM::param_ptr> par_vector = NSM::param_from_r(lP1, lP2, er24_20, k600);
	std::vector<NSM::Pixel> pixes = 
		NSM::dfToPixel(pixdf, light, temperature, pressure, par_vector, dt);

	Rcpp::NumericMatrix result (light.nrow(), light.ncol());
	for(int t = 0; t < light.ncol(); ++t) {
		for(int p = 0; p < pixes.size(); ++p) {
			result(p, t) = pixes.at(p).er(t);
		}
	}
	return result;
}

