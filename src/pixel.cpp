#include "../inst/include/pixel.h"
#include <Rcpp.h>
#include "../inst/include/metabolism.h"
#include "../inst/include/check.h"


/*
	Pixel class
*/

/**
 	METABOLISM FUNCTIONS

 	The following functions compute everything needed for metabolism in a pixel
 	all take current time step as a parameter, and look up internally the relevant information

	simulate_os: Run the one-station simulation for _nt time steps, solving the dDOdt_os
		equation; by default will cache daily gpp and in situ er
	dDOdt_os: Time derivative of dissolved oxygen at time step t for a single station; 
		requires that the simulation be run up to time t-1 before hand
	gpp: Compute gross primary productivity at a given time step g/(m^3 * day)
	er: Compute ecosystem respiration at a given time step g/(m^3 * day)
	rf: Compute reaeration flux at a given time step g/(m^3 * day)
*/
void NSM::Pixel::simulate_os(bool cache) {
	if(cache) {
		_gpp_days = std::vector<double> (ndays(), 0);
		_er_days = std::vector<double> (ndays(), 0);
	}

	for(_timestep = 1; _timestep < _nt; ++_timestep) {
		_DO.at(_timestep) = _DO.at(_timestep - 1) + dDOdt_os(_timestep, cache) * _dt;

	}
}

double NSM::Pixel::dDOdt_os (int t, bool cache) {

	// can't compute for time t=0, requires a previous time step
	if(t < 1)
		throw std::range_error("Time step for dDOdt must be >= 1");

	// can't compute for time past where the simulation has actually run
	if(t > _timestep + 1)
		throw std::runtime_error("Attempted to compute dDOdt at t=" + 
			std::to_string(t) + " but current time step is ct=" + std::to_string(_timestep));

	double cur_gpp = gpp(t-1);
	double cur_er = er(t-1);
	double cur_rf = rf(t-1);

	if(cache) {
		_gpp_days.at(day_from_t(t)) += _dt * cur_gpp / (24 * 60);
		_er_days.at(day_from_t(t)) += _dt * cur_er / (24 * 60);
	}

	// convert everything to be in g/(m^3 * min)
	return ((cur_gpp + cur_er + cur_rf) / _depth) / (24 * 60); 
}

double NSM::Pixel::gpp(int t) {
	return NSM::gpp(_light.at(t), _pars->lP1, _pars->lP2);
}

double NSM::Pixel::er(int t) {
	return NSM::er(_temperature.at(t), _pars->er24_20);
}

double NSM::Pixel::rf(int t) {
	return NSM::reaeration(_temperature.at(t), _pressure.at(t), _DO.at(t), _pars->k600);
}

/**
 	METABOLISM RESULTS FUNCTIONS

 	These public members compute and return the results of the simulation

	daily_gpp: Gross primary productivity, g/(m^3 * day), for each day of the simulation
	daily_insitu_er: Ecosystem respiration, g/(m^3 * day, for each day of the simulation
	do_history: Dissolved oxygen at each time step
*/

const std::vector<double> & NSM::Pixel::daily_gpp() {
	if(_timestep < _nt || _gpp_days.empty())
		throw std::range_error("Run simulation with cache=true before requesting daily gpp");
	return _gpp_days;
}

const std::vector<double> & NSM::Pixel::daily_er() {
	if(_timestep < _nt || _er_days.empty())
		throw std::range_error("Run simulation with cache=true before requesting daily er");
	return _er_days;

}

const std::vector<double> & NSM::Pixel::do_history() {
	if(_timestep < _nt)
		throw std::range_error("Run simulation before requesting the do_history");
	return _DO;
}




NSM::Pixel::Pixel(const param_ptr &pars, double dt, double nt, double depth, double DO_init,
			const std::vector<double> &light, 
			const std::vector<double> &temperature, const std::vector<double> &pressure) : 
		_nt(nt), _dt(dt), _timestep(0), _pars(pars), _depth(depth), _light(light), 
		_temperature(temperature), _pressure(pressure) , _DO(std::vector<double> (nt, -1))
{
	_DO.at(0) = DO_init;

	// input validation needs to happen here
	// still missing: depth, vector length given nt
	check_light(Rcpp::wrap(_light));
	check_DO(Rcpp::wrap(std::vector<double> (_DO.begin(), _DO.begin()+1)));
	check_pressure(Rcpp::wrap(_pressure));

}

/*
	Ways to construct from R objects
*/
std::vector<NSM::Pixel> NSM::dfToPixel(const Rcpp::DataFrame &pixDf, 
		const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature,
		const Rcpp::NumericMatrix &pressure, const std::vector<param_ptr> &pars, 
		double dt) {

	int nt = light.ncol();

	// dimensions checking
	{
		int npix = pixDf.nrow();
		if(light.nrow() != npix)
			throw std::length_error("nrow(light) != nrow(pixdf)");
		if(temperature.nrow() != npix)
			throw std::length_error("nrow(temperature) != nrow(pixdf)");
		if(pressure.nrow() != npix)
			throw std::length_error("nrow(pressure) != nrow(pixdf)");
		if(pars.size() != npix)
			throw std::length_error("length(pars) != nrow(pixdf)");

		if(temperature.ncol() != nt)
			throw std::length_error("ncol(light) != ncol(temperature)");
		if(pressure.ncol() != nt)
			throw std::length_error("ncol(light) != ncol(pressure)");
	}



	std::vector<NSM::Pixel> result;
	Rcpp::NumericVector z = pixDf["depth"];
	Rcpp::NumericVector DO = pixDf["DO_i"];	
	for(int i = 0; i < pixDf.nrow(); ++i) {
		Rcpp::NumericVector lvect = light.row(i);
		std::vector<double> l = Rcpp::as<std::vector<double> >(lvect);
		Rcpp::NumericVector tvect = temperature.row(i);
		std::vector<double> t = Rcpp::as<std::vector<double> >(tvect);
		Rcpp::NumericVector pvect = pressure.row(i);
		std::vector<double> p = Rcpp::as<std::vector<double> >(pvect);
		result.push_back(NSM::Pixel (pars.at(i), dt, nt, z(i), DO(i), l, t, p));
	}
	return(result);
}
