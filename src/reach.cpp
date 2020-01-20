#include "../inst/include/reach.h"
#include "../inst/include/params.h"
#include "../inst/include/pixel.h"
#include <Rcpp.h>

using std::vector;
using std::shared_ptr;

// constructor for the case where light, pressure, and temperature are constant across the reach
// NSM::Reach::Reach(const Rcpp::NumericVector &light, const Rcpp::NumericVector &temperature, 
// 		const Rcpp::NumericVector &pressure, const NSM::Params & pars) : Reach(pars)
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
NSM::Reach::Reach(const Rcpp::DataFrame &pixData, const Rcpp::NumericMatrix &light, 
		const Rcpp::NumericMatrix &temperature, 
		const Rcpp::NumericMatrix &pressure, const NSM::Params & pars, double lateral_do = 0,
		double dt = 10) : Reach(pars)
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
	Rcpp::NumericVector pQ = pixData["discharge"];
	Rcpp::NumericVector pZ = pixData["depth"];
	Rcpp::NumericVector pW = pixData["width"];
	Rcpp::NumericVector pDX = pixData["length"];
	Rcpp::NumericVector pDO = pixData["DO_i"];
	for(int p = 0; p < npix; ++p) {
		const Rcpp::NumericVector lv = light.row(p);
		const Rcpp::NumericVector tv = temperature.row(p);
		const Rcpp::NumericVector pv = pressure.row(p);
		shared_ptr<vector<double> > L (new vector<double>(Rcpp::as<vector<double> >(lv)));
		shared_ptr<vector<double> > T (new vector<double>(Rcpp::as<vector<double> >(tv)));
		shared_ptr<vector<double> > P (new vector<double>(Rcpp::as<vector<double> >(pv)));

		shared_ptr<NSM::Pixel> pix (new NSM::Pixel(_pars, dt, nt, pid.at(p), 
			pQ.at(p), pZ.at(p), pW.at(p), pDX.at(p), pDO.at(p), _lateral_do, L, T, P));
		_pixels[Rcpp::as<int>(pid)] = pix;
	}

	// second pass, set topology
}


NSM::Reach::Reach(const NSM::Params & p) : Reach()
{
	_pars = (shared_ptr<NSM::Params> (new NSM::Params(p)));
}

NSM::Reach::Reach() : _light(new vector<double>), _temperature(new vector<double>),
	_pressure(new vector<double>)
{}
