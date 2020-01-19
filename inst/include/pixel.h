#ifndef pixel_h
#define pixel_h

#include <memory>
#include <vector>
#include <Rcpp.h>
#include "params.h"

namespace NSM {

	class Pixel {
	 	int _nt;  // number of time steps
		int _dt;	// size of time steps, in minutes
		int _timestep;

		param_ptr _pars;
	// 	int id;
	// 	double discharge;
	// 	// double velocity;
		double _depth;
	// 	// double width;
	// 	// double length;

	// 	// eventually provide a friend class LightArray that does the heavy lifting;
	// 	// accessed with light(id, time); thus can abstract all the spatial crap
		std::vector<double> _light; 
		std::vector<double> _temperature; 
		std::vector<double> _pressure;
		std::vector<double> _DO;

		std::vector<double> _gpp_days;
		std::vector<double> _er_days;

	// 	// eventually
	// 	// std::shared_ptr<Pixel> downstream;
	// 	// std::vector<std::shared_ptr<Pixel> > upstream;
		double rf(int t);

	public:
		void simulate_os(bool cache = true);
		double dDOdt_os(int t, bool cache = false);
		double gpp(int t);
		double er(int t);
		int ndays() {return (_nt * _dt) / (24*60);}
		int day_from_t(int t) {return (t * _dt) / (24*60);}

		std::vector<double> daily_gpp();
		std::vector<double> daily_insitu_er();
		std::vector<double> do_history();

		Pixel (const param_ptr &pars, double dt, double nt, double depth, double DO_init, 
			const std::vector<double> &light, 
			const std::vector<double> &temperature, const std::vector<double> &pressure);
	};


	/**
		* Create a vector of pixels from a data.frame and associated data
		* @param pixDf Data.frame giving pixel-specific data
		* @param light Matrix of light values; one row per pixel, one column per time step
		* @param temperature Matrix of water temperature, same dims as light
	  	* @param pressure Matrix of pressure readings, one row per pixel, one column per day
	  	* @param pars vector of pointers to NSM::Params object giving parameters relevant for 
	  			each pixel
	  	* @param dt Size (in minutes) of a time step, defaults to 10 min
	 	* @param nt Number of time steps, defaults to 1 day assuming dt=10

		* @return Vector of pixels
	*/
	std::vector<Pixel> dfToPixel(const Rcpp::DataFrame &pixDf, 
		const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature, 
		const Rcpp::NumericMatrix &pressure, const std::vector<param_ptr> &par, 
		double dt = 10, double nt = (24*60) / 10);

}

#endif
