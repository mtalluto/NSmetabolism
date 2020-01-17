#include <Rcpp.h>
#include <vector>
#include <exception>
#include <string>

using std::range_error;
using std::string;

#include "../inst/include/check.h"

void check_light(const Rcpp::NumericVector &light) {
	check_ge_zero(light, "Light");
}


void check_ge_zero(const Rcpp::NumericVector &dat, string par) {
	for(const auto &d : dat) {
		if(d < 0)
			throw range_error("parameter " + par + " must be >= 0");
	}
}

