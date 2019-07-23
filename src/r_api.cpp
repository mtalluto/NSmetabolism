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
//' @param upDO Upstream dissolved oxygen concentration (mg/L)
//' @param latQ Lateral input discharge (m^3/s)
//' @param latDO Lateral input DO concentration(mg/L)
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

	check_DO(DO);
	check_Q(Q);

	return NSM::oxyMassFlux(DO, Q);
}

//' Compute advective flux, in mg/(m^3 * min)
//' @param inputDOFlux dissolved oxygen mass of input water
//' @param DOconc dissolved oxygen concentration
//' @param data Fixed site data, including at least Q, area, and dx
//' @return transported DO concentration
//' 
// [[Rcpp::export]]
double computeAdvection(double inputDOFlux, double DOconc, const Rcpp::NumericVector &data) {
	check_DO(std::vector<double> {inputDOFlux, DOconc});
	check_site_data(data);

	return NSM::advection();
}