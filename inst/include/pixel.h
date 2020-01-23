#ifndef pixel_h
#define pixel_h

#include <memory>
#include <vector>
#include <unordered_set>
#include <Rcpp.h>

namespace NSM {
	class Params;
	class RData;
	class Timeseries;
	class DailyTimeseries;

	class Pixel {
		size_t _timestep;

		const std::shared_ptr<Params> _pars;
		const std::shared_ptr<RData> _data;
		const std::shared_ptr<Timeseries> _ts_data;
		std::unique_ptr<Timeseries> _DO;
		std::unique_ptr<DailyTimeseries> _daily_totals;

	 	std::unordered_set<std::shared_ptr<Pixel> > _upstream;
	
		double _input_flux(size_t t) const;
		double _ox_mass_flux(size_t t) const;

		double rf(size_t t) const;

	public:
		void add_upstream(std::shared_ptr<Pixel> p) {_upstream.insert(p);};
		void simulate();
		double dDOdt(size_t t, bool cache = false);
		void simulate_os();
		double dDOdt_os(size_t t, bool cache = false);
		double gpp(size_t t) const;
		double er(size_t t) const;
		double advection(size_t t) const;
		unsigned dt() const; // size of time step, minutes
		size_t nt() const; // number of time steps of size dt
		size_t ndays() const; // number of days
		int id() const;
		const std::unordered_set<std::shared_ptr<Pixel> > &upstream() const {return _upstream;}

		std::vector<double> daily_gpp() const;
		std::vector<double> daily_er() const;
		std::vector<double> do_history() const;

		Pixel (std::shared_ptr<Params> pars, std::shared_ptr<RData> data, 
			const Timeseries &ts);
		Pixel (std::shared_ptr<Params> pars, std::shared_ptr<RData> data, 
			std::shared_ptr<Timeseries> ts);
	};


	/**
		* Create a vector of pixels from a data.frame and associated data
		* @param pixDf Data.frame giving pixel-specific data
		* @param light Matrix of light values; one row per pixel, one column per time step
		* @param temperature Matrix of water temperature, same dims as light
	  	* @param pressure Matrix of pressure readings, one row per pixel, one column per day
	  	* @param pars vector of pointers to NSM::Params object giving parameters relevant for 
	  			each pixel
	  	* @param dt Size (in minutes) of a time step, defaults to 10 min; number of time steps
	 		will be inferred from the number of columns

		* @return Vector of pixels
	*/
	std::vector<Pixel> dfToPixel(const Rcpp::DataFrame &pixDf, 
		const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature, 
		const Rcpp::NumericMatrix &pressure, const std::vector<std::shared_ptr<Params> > &par, 
		double dt = 10);

}

#endif
