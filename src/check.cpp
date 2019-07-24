#include <Rcpp.h>
#include <vector>
#include <string>
#include <exception>
#include "../inst/include/check.h"

using std::vector;
using std::range_error;
using std::string;

void check_DO(const Rcpp::NumericVector &DO) {
	check_ge_zero(DO, "DO concentration");
}

void check_Q(const Rcpp::NumericVector &Q) {
	check_ge_zero(Q, "Discharge");
}

void check_pressure(const Rcpp::NumericVector &P) {
	check_ge_zero(P, "Pressure");
}

void check_site_data(const Rcpp::NumericVector &site) {
	check_ge_zero(Rcpp::as<Rcpp::NumericVector>(site["Q"]), "Q");
	Rcpp::CharacterVector dims ({"area", "dx"});
	Rcpp::NumericVector dimDat (Rcpp::as<Rcpp::NumericVector>(site[dims]));
	check_positive(dimDat, "Dimensions (area, length, depth)");
}

void check_light(const Rcpp::NumericVector &light) {
	check_ge_zero(light, "PAR");
}

void check_k(const Rcpp::NumericVector &k) {
	check_ge_zero(k, "k");
}


void check_ge_zero(const Rcpp::NumericVector &dat, string par) {
	for(const auto &d : dat) {
		if(d < 0)
			throw range_error("parameter " + par + " must be >= 0");
	}
}

void check_positive(const Rcpp::NumericVector &dat, string par) {
	for(const auto &d : dat) {
		if(d <= 0)
			throw range_error("parameter " + par + " must be positive");
	}
}
