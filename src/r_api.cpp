#include <Rcpp.h>
#include <vector>
#include "../inst/include/r_api.h"
#include "../inst/include/metabolism.h"
#include "../inst/include/check.h"

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


//' Compute gross primary productivity from light (g O2 / [m^2 * day]) using pixels
//' @param pixdf The pixel data frame to construct the pixels
//' @param light Matrix of light readings; pixels in rows, light readings in columns
//' @param lP1 The log of the slope of the PI curve
//' @param lP2 The log of the saturation term of the PI curve
//' @return Matrix; computed GPP for each pixel at each time step
Rcpp::NumericMatrix pixelGPP(const Rcpp::DataFrame &pixdf, const Rcpp::NumericMatrix &light, 
		double lP1, double lP2) {
	for(i = 0; i < light.nrow(); ++i)
		check_light(light.row(i));
	if(light.nrow() != pixdf.nrow())
		throw std::length_error("pixdf and light must have same number of rows");

	Rcpp::NumericMatrix result (light.nrow(), light.ncol());
	std::vector<NSM::Pixel> pixes;
	for(i = 0; i < pixdf.nrow(); ++i) {
		pixes.push_back(Pixel(..., light.row(i)));
	}

	for(int t = 0; t < light.ncol(); ++t) {
		for(p in pixes.size())
		result(p, t) = pixes.at(p).gpp(t);
	}
}
