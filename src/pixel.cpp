#include "../inst/include/pixel.h"
#include "../inst/include/params.h"
#include "../inst/include/metabolism.h"
#include "../inst/include/check.h"
#include <Rcpp.h>

using std::shared_ptr;
using std::vector;
using std::make_shared;

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
void NSM::Pixel::simulate() {
	for(_timestep = 1; _timestep < nt(); ++_timestep) {
		_DO->at("DO", _timestep) = _DO->at("DO", _timestep - 1) + dDOdt(_timestep, true) * dt();
	}
}

double NSM::Pixel::dDOdt (size_t t, bool cache) {
	return advection(t) + dDOdt_os(t, cache); 
}


void NSM::Pixel::simulate_os() {
	for(_timestep = 1; _timestep < nt(); ++_timestep) {
		_DO->at("DO", _timestep) = _DO->at("DO", _timestep - 1) + dDOdt_os(_timestep, true) * dt();
	}
}

double NSM::Pixel::dDOdt_os (size_t t, bool cache) {
	// can't compute for time past where the simulation has actually run
	if(t > _timestep + 1)
		throw std::runtime_error("Attempted to compute dDOdt at t=" + 
			std::to_string(t) + " but current time step is ct=" + std::to_string(_timestep));

	double cur_gpp = gpp(t-1);
	double cur_er = er(t-1);
	double cur_rf = rf(t-1);

	if(cache) {
		_daily_totals->at_min("gpp", t * dt()) += dt() * cur_gpp / NSM::MIN_PER_DAY;
		_daily_totals->at_min("er", t * dt()) += dt() * cur_er / NSM::MIN_PER_DAY;
	}

	// convert everything to be in g/(m^3 * min)
	return ((cur_gpp + cur_er + cur_rf) / _data->depth) / NSM::MIN_PER_DAY; 
}

double NSM::Pixel::gpp(size_t t) const {
	return NSM::gpp(_ts_data->at("light", t), _pars->lP1, _pars->lP2);
}

double NSM::Pixel::er(size_t t) const {
	return NSM::er(_ts_data->at("temperature", t), _pars->er24_20);
}

double NSM::Pixel::rf(size_t t) const {
	return NSM::reaeration(_ts_data->at("temperature", t), _ts_data->at("pressure", t), 
		_DO->at("DO", t), _pars->k600);
}


/**
	* Compute advective concentration flux (g/[m^3 * min])
	* @return double; dissolved oxygen concentration flux (g/[m^3 * min])
*/

double NSM::Pixel::advection(size_t t) const {

	// Input DO mass flux (g/min) from upstream and lateral input for previous time step
	double inputFlux = _input_flux(t - 1);

	// Output DO mass flux (g/min) for the previous time step 
	double outputFlux = _ox_mass_flux(t - 1);

	return (-1/_data->cs_area()) * (outputFlux - inputFlux)/_data->length;
}


/**
	* Compute oxygen flux from the pixel at a particular time (g/min)
*/
double NSM::Pixel::_ox_mass_flux(size_t t) const {
	return NSM::oxygen_flux(_DO->at("DO", t), _data->discharge);
}


/**
	* Compute input flux (g/min) by summing across all sources
 	* @return total input flux (g/min)
*/
double NSM::Pixel::_input_flux(size_t t) const {
	double inFlux = 0;
	double inQ = 0;
	if(!_upstream.empty()) {
		for(const auto &us : _upstream) {
			inFlux += us->_ox_mass_flux(t);
			inQ += us->_data->discharge;
		}
	}

	double lateralQ = _data->discharge - inQ;
	if(lateralQ < 0)
		lateralQ = 0;

	inFlux += _data->lateral_do * lateralQ;
	return inFlux;
}


/**
 	METABOLISM RESULTS FUNCTIONS

 	These public members compute and return the results of the simulation

	daily_gpp: Gross primary productivity, g/(m^3 * day), for each day of the simulation
	daily_insitu_er: Ecosystem respiration, g/(m^3 * day, for each day of the simulation
	do_history: Dissolved oxygen at each time step
*/

vector<double> NSM::Pixel::daily_gpp() const {
	if(_timestep < nt())
		throw std::out_of_range("Run simulation with cache=true before requesting daily gpp");
	return _daily_totals->at("gpp");
}

vector<double> NSM::Pixel::daily_er() const {
	if(_timestep < nt())
		throw std::out_of_range("Run simulation with cache=true before requesting daily er");
	return _daily_totals->at("er");
}

vector<double> NSM::Pixel::do_history() const {
	if(_timestep < nt())
		throw std::out_of_range("Run simulation before requesting the do_history");
	return _DO->at("DO");
}



/**
 	HELPER FUNCTIONS

 	dt: size of time step, in minutes
 	nt: number of time steps
*/
unsigned NSM::Pixel::dt() const {return _ts_data->dt();}
size_t NSM::Pixel::nt() const {return _ts_data->nt();}
size_t NSM::Pixel::ndays() const {return _daily_totals->nt();}
int NSM::Pixel::id() const {return _data->id;}

/*
 	CONSTRUCTORS
*/

/*
	Constructor with a copy of a time series (so individual to this pixel)
*/
NSM::Pixel::Pixel(shared_ptr<Params> pars, shared_ptr<RData> data,
	const NSM::Timeseries &ts) : Pixel(pars, data, make_shared<NSM::Timeseries>(ts))
{ }

/*
	Constructor with shared time series
*/
NSM::Pixel::Pixel (shared_ptr<Params> pars, shared_ptr<RData> data, 
	shared_ptr<NSM::Timeseries> ts) :
		_timestep(0), _pars(pars), _data(data), _ts_data(ts)
{
	_DO = std::make_unique<NSM::Timeseries> 
		(_ts_data->dt(), "DO", vector<double> (_ts_data->nt(), -1));
	_DO->at("DO", 0) = _data->initial_do;

	int nday = (nt() * dt()) / NSM::MIN_PER_DAY;
	if(((nt() * dt()) % NSM::MIN_PER_DAY) != 0) {
		nday++;
		Rcpp::Rcerr << "warning: number of minutes to run simulation is not an even " << 
			"number of days; gpp and er daily totals will be inaccurate. To prevent this " <<
			"warning, be sure to start the model at midnight of the first day and end at " <<
			"midnight - dt of the last day; or use an even number of 24-h periods.\n";
	}
	_daily_totals = std::make_unique<NSM::DailyTimeseries>();
	_daily_totals->insert("gpp", vector<double> (nday, 0));
	_daily_totals->insert("er", vector<double> (nday, 0));

	// input validation needs to happen here
	// items remaining to check:
	//    dt, nt, discharge, depth, width, length, lateral, 
	check_light(Rcpp::wrap(_ts_data->at("light")));
	check_DO(Rcpp::wrap(vector<double> (1, _data->initial_do)));
	check_pressure(Rcpp::wrap(_ts_data->at("pressure")));
}



/*
	Ways to construct from R objects
*/
vector<NSM::Pixel> NSM::dfToPixel(const Rcpp::DataFrame &pixDf, 
		const Rcpp::NumericMatrix &light, const Rcpp::NumericMatrix &temperature,
		const Rcpp::NumericMatrix &pressure, 
		const vector<shared_ptr<NSM::Params> > &pars, 
		double dt) {

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
	}

	vector<NSM::Pixel> result;
	Rcpp::NumericVector z = pixDf["depth"];
	Rcpp::NumericVector DO = pixDf["DO_i"];	
	for(int i = 0; i < pixDf.nrow(); ++i) {
		auto ts = std::make_shared<NSM::Timeseries> (dt, "light", light.row(i));
		ts->insert("temperature", temperature.row(i));
		ts->insert("pressure", pressure.row(i));
		auto dat = make_shared<RData> (z(i), DO(i));
		result.push_back(NSM::Pixel (pars.at(i), dat, ts));
	}
	return(result);
}
