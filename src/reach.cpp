#include "../inst/include/reach.h"
#include "../inst/include/params.h"
#include "../inst/include/pixel.h"
#include <Rcpp.h>

using std::vector;
using std::shared_ptr;
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
		double lateral_do = 0, double dt = 10) : Reach(pars)
{
	_lateral_do = lateral_do;

	int npix = light.nrow();
	int nt = light.ncol();

	if(temperature.nrow() != npix || pressure.nrow() != npix)
		throw std::length_error("Light, temperature, and pressure must have the same nrow");
	if(temperature.ncol() != nt || pressure.ncol() != nt)
		throw std::length_error("Light, temperature, and pressure must have the same ncol");

	// First pass, set up pixel data frame
	Rcpp::IntegerVector pid = pixData["id"];
	NumericVector pQ = pixData["discharge"];
	NumericVector pZ = pixData["depth"];
	NumericVector pW = pixData["width"];
	NumericVector pDX = pixData["length"];
	NumericVector pDO = pixData["DO_i"];
	for(int p = 0; p < npix; ++p) {
		const NumericVector lv = light.row(p);
		const NumericVector tv = temperature.row(p);
		const NumericVector pv = pressure.row(p);
		shared_ptr<vector<double> > L (new vector<double>(Rcpp::as<vector<double> >(lv)));
		shared_ptr<vector<double> > T (new vector<double>(Rcpp::as<vector<double> >(tv)));
		shared_ptr<vector<double> > P (new vector<double>(Rcpp::as<vector<double> >(pv)));

		shared_ptr<NSM::Pixel> pix (new NSM::Pixel(_pars, dt, nt, 
			pQ.at(p), pZ.at(p), pW.at(p), pDX.at(p), pDO.at(p), _lateral_do, L, T, P));
		_pixels[pid.at(p)] = pix;
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


NSM::Reach::Reach(const NSM::Params & p) : Reach()
{
	_pars = (shared_ptr<NSM::Params> (new NSM::Params(p)));
}

NSM::Reach::Reach() : _light(new vector<double>), _temperature(new vector<double>),
	_pressure(new vector<double>)
{}
