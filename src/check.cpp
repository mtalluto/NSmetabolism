#include <Rcpp.h>
#include <vector>
#include <exception>
#include "../inst/include/check.h"

using std::vector;
using std::range_error;

void check_DO(const vector<double> &DO) {
	for(const auto & iDO : DO) {
		if(iDO < 0)
			throw range_error("DO concentrations must be >= 0");
	}
}

void check_Q(const vector<double> &Q) {
	for(const auto & iQ : Q) {
		if(iQ < 0)
			throw range_error("Q must be >= 0");
	}	
}


check_site_data(Rcpp::NumericVector site) {
	if(site["Q"] < 0)
		throw std::range_error("Q must be >= 0");

	if(site["area"] <= 0)
		throw std::range_error("Area must positive");

	if(site["dx"] <= 0)
		throw std::range_error("dx must positive");

}
