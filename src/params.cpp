#include "../inst/include/params.h"
#include <Rcpp.h>

using std::shared_ptr;

NSM::Params::Params(const NSM::Params &p) : lP1(p.lP1), lP2(p.lP2), er24_20(p.er24_20), 
	k600(p.k600)
{}

NSM::Params::Params(double lp1, double lp2, double er, double k) : lP1(lp1), lP2(lp2), 
	er24_20(er), k600(k)
{}

// NSM::Params& NSM::Params::operator= (const Params& rhs) {
// 	lP1 = rhs.lP1;
// 	lP2 = rhs.lP2;
// 	er24_20 = rhs.er24_20;
// 	k600 = rhs.k600;
// 	return *this;
// }


std::vector<shared_ptr<NSM::Params> > NSM::param_from_r(const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
		const Rcpp::NumericVector &k600)
{
	int ni = lP1.size();
	if(lP2.size() != ni || er24_20.size() != ni || k600.size() != ni)
		throw std::length_error("All parameters must have the same length");

	std::vector<shared_ptr<NSM::Params> > par_vector;
	for(int i = 0; i < ni; ++i)
		par_vector.push_back(shared_ptr<NSM::Params> 
			(new NSM::Params {lP1.at(i), lP2.at(i), er24_20.at(i), k600.at(i)}));
	return par_vector;
}

