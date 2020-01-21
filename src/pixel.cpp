#include "../inst/include/pixel.h"
#include "../inst/include/params.h"
#include "../inst/include/metabolism.h"
#include "../inst/include/check.h"
#include <Rcpp.h>

/*
	Pixel class
*/

/**
 	METABOLISM FUNCTIONS

 	The following functions compute everything needed for metabolism in a pixel
 	all take current time step as a parameter, and look up internally the relevant information

	simulate_os: Run the n-station simulation for _nt time steps, solving the dDOdt
		equation; by default will cache daily gpp and in situ er
	dDOdt: Time derivative of dissolved oxygen at time step t considering upstream/lateral input; 
		requires that the simulation be run up to time t-1 before hand
	simulate_os: Run the one-station simulation for _nt time steps, solving the dDOdt_os
		equation; by default will cache daily gpp and in situ er
	dDOdt_os: Time derivative of dissolved oxygen at time step t for a single station; 
		requires that the simulation be run up to time t-1 before hand
	gpp: Compute gross primary productivity at a given time step g/(m^3 * day)
	er: Compute ecosystem respiration at a given time step g/(m^3 * day)
	rf: Compute reaeration flux at a given time step g/(m^3 * day)
*/
void NSM::Pixel::simulate(bool cache) {
	if(cache) {
		_gpp_days = std::vector<double> (ndays(), 0);
		_er_days = std::vector<double> (ndays(), 0);
	}

	for(_timestep = 1; _timestep < _nt; ++_timestep) {
		_DO.at(_timestep) = _DO.at(_timestep - 1) + dDOdt(_timestep, cache) * _dt;
	}
}

double NSM::Pixel::dDOdt (int t, bool cache) {
	return advection(t) + dDOdt_os(t, cache); 
}


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
	return NSM::gpp(_light->at(t), _pars->lP1, _pars->lP2);
}

double NSM::Pixel::er(int t) {
	return NSM::er(_temperature->at(t), _pars->er24_20);
}

double NSM::Pixel::rf(int t) {
	return NSM::reaeration(_temperature->at(t), _pressure->at(t), _DO.at(t), _pars->k600);
}


/**
	* Compute advective concentration flux (g/[m^3 * min])
	* @return double; dissolved oxygen concentration flux (g/[m^3 * min])
*/

double NSM::Pixel::advection(int t) {

	// Input DO mass flux (g/min) from upstream and lateral input for previous time step
	double inputFlux = _input_flux(t - 1);

	// Output DO mass flux (g/min) for the previous time step 
	double outputFlux = _ox_mass_flux(t - 1);

	return (-1/_area()) * (outputFlux - inputFlux)/_length;
}


/**
	* Compute oxygen flux from the pixel at a particular time (g/min)
*/
double NSM::Pixel::_ox_mass_flux(int t) const {
	return NSM::oxygen_flux(_DO.at(t), _discharge);
}


/**
	* Compute input flux (g/min) by summing across all sources
 	* @return total input flux (g/min)
*/
double NSM::Pixel::_input_flux(int t) const {
	double inFlux = 0;
	double inQ = 0;
	if(!_upstream.empty()) {
		for(const auto &us : _upstream) {
			inFlux += us->_ox_mass_flux(t);
			inQ += us->_discharge;
		}
	}

	double lateralQ = _discharge - inQ;
	if(lateralQ < 0)
		lateralQ = 0;

	inFlux += _lateral_do * lateralQ;
	return inFlux;
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





/*
 	CONSTRUCTORS
*/

/*
	Simple constructor for a single pixel, meant for one-station model
*/
NSM::Pixel::Pixel(const std::shared_ptr<Params> &pars, double dt, double nt, double depth, 
		double DO_init, const std::vector<double> &light, const std::vector<double> &temperature,
		const std::vector<double> &pressure) : 
	Pixel(pars, dt, nt, 0, depth, 0, 0, DO_init, 0, 
		std::make_shared<std::vector<double> >(light), 
		std::make_shared<std::vector<double> >(temperature), 
		std::make_shared<std::vector<double> >(pressure))
{ }

/*
	Constructor for the n-station model, which has additional data requirements
*/
NSM::Pixel::Pixel (const std::shared_ptr<Params> &pars, double dt, double nt, 
			double discharge, double depth, double width, double length, double DO_init, 
			double lateral, std::shared_ptr<std::vector<double> > light, 
			std::shared_ptr<std::vector<double> > temperature, 
			std::shared_ptr<std::vector<double> > pressure) :
		_nt(nt), _dt(dt), _timestep(0), _pars(pars), _discharge(discharge), _depth(depth),
		_width(width), _length(length), _lateral_do(lateral), _light(light),
		_temperature(temperature), _pressure(pressure) , _DO(std::vector<double> (nt, -1))
{
	_DO.at(0) = DO_init;

	// input validation needs to happen here
	// items remaining to check:
	//    dt, nt, discharge, depth, width, length, lateral, 
	check_light(Rcpp::wrap(*_light));
	check_DO(Rcpp::wrap(std::vector<double> (_DO.begin(), _DO.begin()+1)));
	check_pressure(Rcpp::wrap(*_pressure));
}



/*
	Ways to construct from R objects
*/
std::vector<NSM::Pixel> NSM::dfToPixel(const Rcpp::DataFrame &pixDf, 
		const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature,
		const Rcpp::NumericMatrix &pressure, 
		const std::vector<std::shared_ptr<NSM::Params> > &pars, 
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
