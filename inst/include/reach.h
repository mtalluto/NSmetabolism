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
		std::shared_ptr<std::vector<double> > _light;
		std::shared_ptr<std::vector<double> > _temperature;
		std::shared_ptr<std::vector<double> > _pressure;
		std::shared_ptr<Params> _pars;
		double _lateral_do;

		// map with pixelID as key
		std::map<int, std::shared_ptr<Pixel> > _pixels;

		// pointers to the pixels at the top and bottom of the reach
		std::shared_ptr<Pixel> _bottom;
		std::shared_ptr<Pixel> _top;


		Reach();
		Reach(const Params & p);
	public:
		// Reach(const Rcpp::NumericVector &light, const Rcpp::NumericVector &temperature, 
		// 	const Rcpp::NumericVector &pressure, const Params &pars);
		Reach(const Rcpp::DataFrame &pixData, const Rcpp::NumericMatrix &light, 
			const Rcpp::NumericMatrix &temperature, const Rcpp::NumericMatrix &pressure, 
			const NSM::Params & pars, const Rcpp::NumericMatrix &topology, int bottom, 
			int top, double lateral_do, double dt);
	};
}

#endif
