#ifndef params_h
#define params_h

#include <memory>

/**
 * @brief      Class for holding parameters for the metabolism model 
*/
namespace NSM {
	class Params {
	public:
		double lP1;
		double lP2;
		double er24_20;
		double k600;
	};

	typedef std::shared_ptr<Params> param_ptr;
	std::vector<param_ptr> param_from_r(const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
		const Rcpp::NumericVector &k600);
}

#endif
