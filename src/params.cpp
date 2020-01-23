#include "../inst/include/params.h"
#include <Rcpp.h>

using std::shared_ptr;
using Rcpp::NumericVector;

NSM::Params::Params(const NSM::Params &p) : lP1(p.lP1), lP2(p.lP2), er24_20(p.er24_20), 
	k600(p.k600)
{}

NSM::Params::Params(double lp1, double lp2, double er, double k) : lP1(lp1), lP2(lp2), 
	er24_20(er), k600(k)
{}

NSM::Params::Params(const NumericVector &p) : Params(p("lP1"), p("lP2"), p("er24_20"), p("k600"))
{}

// NSM::Params& NSM::Params::operator= (const Params& rhs) {
// 	lP1 = rhs.lP1;
// 	lP2 = rhs.lP2;
// 	er24_20 = rhs.er24_20;
// 	k600 = rhs.k600;
// 	return *this;
// }

const std::set<std::string> NSM::Timeseries::_allowed_variables(
	{"light", "temperature", "pressure", "DO", "gpp", "er"});

void NSM::Timeseries::insert(std::string var, std::vector<double> data) {
	if(_allowed_variables.find(var) == _allowed_variables.end())
		throw std::out_of_range("Tried to insert unknown timeseries variable: " + var);
	if(_data.find(var) != _data.end())
		throw std::out_of_range("Tried to insert an already existing Timeseries variable: " + var);
	if(!empty() && data.size() != nt())
		throw std::out_of_range("Tried to insert a time series of size " + 
			std::to_string(data.size()) + " into a timeseries with a size of " + 
			std::to_string(nt()));
	_data.insert(std::pair<std::string, std::vector<double> > (var, data));
}

void NSM::Timeseries::insert(std::string var, const Rcpp::NumericVector &data) {
	insert(var, Rcpp::as<std::vector<double> >(data));
}

double NSM::Timeseries::at(std::string var, size_t time) const {
	return _data.at(var).at(time);
}

std::vector<double> NSM::Timeseries::at(std::string var) const {
	return _data.at(var);
}

double& NSM::Timeseries::at(std::string var, size_t time) {
	return _data.at(var).at(time);
}


// same as at, but accepts the time elapsed in minutes, rather than in days
double& NSM::DailyTimeseries::at_min(std::string var, size_t time_min) {
	return _data.at(var).at(time_min / NSM::MIN_PER_DAY);
}


NSM::Timeseries::Timeseries(unsigned dt) : _dt(dt)
{}

NSM::Timeseries::Timeseries(unsigned dt, std::string var, std::vector<double> data) : 
NSM::Timeseries(dt)
{
	insert(var, data);
}

NSM::Timeseries::Timeseries(unsigned dt, std::string var, const Rcpp::NumericVector &data) :
NSM::Timeseries(dt, var, Rcpp::as<std::vector<double> >(data))
{}




std::vector<shared_ptr<NSM::Params> > NSM::param_from_r(const Rcpp::NumericVector &lP1, 
		const Rcpp::NumericVector &lP2, const Rcpp::NumericVector &er24_20, 
		const Rcpp::NumericVector &k600)
{
	int ni = lP1.size();
	if(lP2.size() != ni || er24_20.size() != ni || k600.size() != ni)
		throw std::length_error("All parameters must have the same length");

	std::vector<shared_ptr<NSM::Params> > par_vector;
	for(int i = 0; i < ni; ++i)
		par_vector.push_back(shared_ptr<NSM::Params> 
			(new NSM::Params {lP1.at(i), lP2.at(i), er24_20.at(i), k600.at(i)}));
	return par_vector;
}


std::vector<shared_ptr<NSM::RData> > NSM::data_from_r(const Rcpp::DataFrame &df) {
	int n = df.nrow();
	Rcpp::IntegerVector pid = df["id"];
	NumericVector pQ = df["discharge"];
	NumericVector pZ = df["depth"];
	NumericVector pW = df["width"];
	NumericVector pDX = df["length"];
	NumericVector lDO = df["lateral"];
	NumericVector pDO = df["DO_i"];
	std::vector<shared_ptr<NSM::RData> > result;
	for(int i = 0; i < n; ++i) {
		auto pDat = std::make_shared<NSM::RData> (pZ.at(i), pDO.at(i), pQ.at(i), pW.at(i), 
			pDX.at(i), lDO.at(i), pid.at(i));
		result.push_back(pDat);
	}
	return result;
}



