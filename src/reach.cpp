#include "../inst/include/reach.h"
#include "../inst/include/params.h"
#include "../inst/include/pixel.h"
#include <Rcpp.h>

using std::vector;
using std::shared_ptr;
using std::make_shared;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// constructor for the case where light, pressure, and temperature are constant across the reach
// NSM::Reach::Reach(const NumericVector &light, const NumericVector &temperature, 
// 		const NumericVector &pressure, const NSM::Params & pars) : Reach(pars)
// {
// 	int nt = light.size();
// 	if(temperature.size() != nt || pressure.size() != nt)
// 		throw std::length_error("Light, temperature, and pressure must have the same length");

// 	for(int i = 0; i < nt; ++i) {
// 		_light->push_back(light.at(i));
// 		_temperature->push_back(temperature.at(i));
// 		_pressure->push_back(pressure.at(i));
// 	}
// }


// constructor for the case where light, pressure, and temperature vary by pixel
NSM::Reach::Reach(const Rcpp::DataFrame &pixData, const NumericMatrix &light, 
		const NumericMatrix &temperature, 
		const NumericMatrix &pressure, const NSM::Params & pars, 
		const NumericMatrix &topology, int bottom, int top, 
		double dt = 10) : Reach(pars)
{
	int npix = light.nrow();
	int nt = light.ncol();

	if(temperature.nrow() != npix || pressure.nrow() != npix)
		throw std::length_error("Light, temperature, and pressure must have the same nrow");
	if(temperature.ncol() != nt || pressure.ncol() != nt)
		throw std::length_error("Light, temperature, and pressure must have the same ncol");

	// First pass, set up pixels
	auto data = data_from_r(pixData);
	for(int p = 0; p < npix; ++p) {
		auto ts = make_shared<NSM::Timeseries> (dt, "light", light.row(p));
		ts->insert("temperature", temperature.row(p));
		ts->insert("pressure", pressure.row(p));

		_pixels[data.at(p)->id] = make_shared<NSM::Pixel> (_pars, data.at(p), ts);
	}

	// second pass, set topology
	_bottom = _pixels.at(bottom);
	_top = _pixels.at(top);
	for(int i = 0; i < topology.nrow(); ++i) {
		int us = topology(i, 0);
		int ds = topology(i, 1);
		_pixels.at(ds)->add_upstream(_pixels.at(us));
	}
}


NSM::Reach::Reach(const NSM::Params & p)
{
	_pars = (shared_ptr<NSM::Params> (new NSM::Params(p)));
}

