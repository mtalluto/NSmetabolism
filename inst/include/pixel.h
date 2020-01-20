#ifndef pixel_h
#define pixel_h

#include <memory>
#include <vector>
#include <Rcpp.h>

namespace NSM {
	class Params;

	class Pixel {
	 	int _nt;  // number of time steps
		int _dt;	// size of time steps, in minutes
		int _timestep;

		std::shared_ptr<Params> _pars;
	// 	int id;
		double _discharge; // m^3/sec
		double _depth; // in m
		double _width; // in m
		double _length; // in m
		double _lateral_do; // in g/m^3

	// 	// eventually provide a friend class LightArray that does the heavy lifting;
	// 	// accessed with light(id, time); thus can abstract all the spatial crap
		std::vector<double> _light; 
		std::vector<double> _temperature; 
		std::vector<double> _pressure;
		std::vector<double> _DO;

		std::vector<double> _gpp_days;
		std::vector<double> _er_days;

	// 	// std::shared_ptr<Pixel> downstream;
	 	std::vector<std::shared_ptr<Pixel> > _upstream;
	
		double _input_flux(int t) const;
		double _ox_mass_flux(int t) const;
		double _area() const {return _depth * _width;}

		double rf(int t);

	public:
		void simulate(bool cache = true);
		double dDOdt(int t, bool cache = false);
		void simulate_os(bool cache = true);
		double dDOdt_os(int t, bool cache = false);
		double gpp(int t);
		double er(int t);
		double advection(int t);
		int ndays() {return (_nt * _dt) / (24*60);}
		int day_from_t(int t) {return (t * _dt) / (24*60);}

		const std::vector<double> & daily_gpp();
		const std::vector<double> & daily_er();
		const std::vector<double> & do_history();

		Pixel (const std::shared_ptr<Params> &pars, double dt, double nt, double depth, 
			double DO_init, const std::vector<double> &light, 
			const std::vector<double> &temperature, const std::vector<double> &pressure);

		// STOP - need to implement this guy, which means changing how Pixel is implemented
		// behind the scenes, needs pointers for light temperature etc
		// don't have to change the ctor above, can just copy the vectors to not break earlier
		// functions
		Pixel (const std::shared_ptr<Params> &pars, double dt, double nt, int id, 
			double discharge, double depth, double width, double length, double DO_init, 
			double lateral, std::shared_ptr<std::vector<double> > light, 
			std::shared_ptr<std::vector<double> > temperature, 
			std::shared_ptr<std::vector<double> > pressure);
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
