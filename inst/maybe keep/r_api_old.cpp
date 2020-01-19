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









