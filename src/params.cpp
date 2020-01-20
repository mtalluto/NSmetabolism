#include "../inst/include/pixel.h"
#include <Rcpp.h>

std::vector<NSM::param_ptr> NSM::param_from_r(const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
		const Rcpp::NumericVector &k600)
{
	int ni = lP1.size();
	if(lP2.size() != ni || er24_20.size() != ni || k600.size() != ni)
		throw std::length_error("All parameters must have the same length");

	std::vector<NSM::param_ptr> par_vector;
	for(int i = 0; i < ni; ++i)
		par_vector.push_back(NSM::param_ptr 
			(new NSM::Params {lP1.at(i), lP2.at(i), er24_20.at(i), k600.at(i)}));
	return par_vector;
}
