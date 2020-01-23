#ifndef reach_h
#define reach_h

#include <unordered_map>
#include <vector>
#include <memory>
#include <Rcpp.h>


namespace NSM {
	class Params;
	class Pixel;
	class Timeseries;

	class Reach {
		// std::shared_ptr<Timeseries> _ts_data;
		std::shared_ptr<Params> _pars;
		bool _topology_built;

		// map with pixelID as key
		std::unordered_map<int, std::shared_ptr<Pixel> > _pixels;

		// pointers to the pixels at the top and bottom of the reach
		std::shared_ptr<Pixel> _bottom;
		std::shared_ptr<Pixel> _top;

		std::unordered_set<std::shared_ptr<Reach> > _upstream;

		void set_topology(const Rcpp::IntegerMatrix &T);


		Reach(std::shared_ptr<Params> p);
	public:
		std::shared_ptr<Pixel> bottom();
		// Reach(const Rcpp::NumericVector &light, const Rcpp::NumericVector &temperature, 
		// 	const Rcpp::NumericVector &pressure, const Params &pars);
		Reach(std::shared_ptr<Params> pars, const Rcpp::DataFrame &pixData, 
			std::vector<std::shared_ptr<Timeseries> > ts_data, 
			const Rcpp::IntegerMatrix &topology);
	};
}

#endif
