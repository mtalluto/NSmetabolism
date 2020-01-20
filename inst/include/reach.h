#ifndef reach_h
#define reach_h

#include <map>
#include <vector>
#include <memory>
#include <Rcpp.h>


namespace NSM {
	class Params;
	class Pixel;

	class Reach {
		// we assume light, water temperature, pressure. 
		// parameters are not variable at the reach scale
		std::shared_ptr<std::vector<double> > _light;
		std::shared_ptr<std::vector<double> > _temperature;
		std::shared_ptr<std::vector<double> > _pressure;
		std::shared_ptr<Params> _pars;
		double _lateral_do;

		// map with pixelID as key
		std::map<int, std::shared_ptr<Pixel> > _pixels;

		Reach();
		Reach(const Params & p);
	public:
		// Reach(const Rcpp::NumericVector &light, const Rcpp::NumericVector &temperature, 
		// 	const Rcpp::NumericVector &pressure, const Params &pars);
		Reach(const Rcpp::DataFrame &pixData, const Rcpp::NumericMatrix &light, 
			const Rcpp::NumericMatrix &temperature, const Rcpp::NumericMatrix &pressure, 
			const NSM::Params & pars, double lateral_do, double dt);
		// Reach(, 
		// const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature, 
		// const Rcpp::NumericMatrix &pressure, const std::vector<std::shared_ptr<Params> > &par, 
		// double dt = 10);
	};
}

#endif
